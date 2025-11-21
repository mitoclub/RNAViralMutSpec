#!/usr/bin/env python3
"""
this script must work but it's old version with some shit
"""
import os
import re
from multiprocessing import Process
from time import sleep

import datetime
from selenium.webdriver.firefox.options import Options as FirefoxOptions
from selenium.webdriver.chrome.options import Options as ChromeOptions
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By
from selenium import webdriver

from base import Parser

GISAID_URL = "https://www.gisaid.org/"
ACCESS_FILE_PATH = "./access_file.txt"
LAST_IDX_PATH = "./last.id"

START_IDX = "EPI_ISL_1044108"
LAST_IDX =  "EPI_ISL_5144811"

SEQ_AMOUNT = 2000
INSERT_CHANK_SIZE=2000


def read_last_idx():
    """ idx that we can start from """
    if not os.path.exists(LAST_IDX_PATH):
        return START_IDX

    with open(LAST_IDX_PATH) as fin:
        last_idx = fin.read()
    return last_idx


def write_last_idx(last_idx):
    assert isinstance(last_idx, str), "ID must be string"
    with open(LAST_IDX_PATH, "w") as fout:
        fout.write(last_idx)


class GisaidParser(Parser):
    def __init__(self, headless=False, wdriver="chrome"):
        super().__init__(headless=headless)
        
        self.read_access_data()
        self.cur_date = "2020-04-28"  # this day downloaded but we will parse next days

    def read_access_data(self, path=ACCESS_FILE_PATH):
        with open(path) as fin:
            login, password = fin.readline().strip().split()
        self.__login = login
        self.__password = password

    def get_total_number(self, driver=None):
        driver = driver or self.driver
        wait = WebDriverWait(driver, 10)
        wait.until(EC.text_to_be_present_in_element(
            (By.CLASS_NAME, 'sys-datatable-info-left'), 'Total'))  # EC.func() return bool

        total_raw = self.find_elem(".sys-datatable-info-left").text
        total = re.search('Total: (.+) viruses', total_raw).groups()[0]
        total = total.replace(',', '')
        return int(total)

    def login_and_open_table(self, url=GISAID_URL):
        self.driver.get(url)
        self.find_and_click(".Login")

        self.find_and_fill("#elogin", self.__login)
        self.find_and_fill("#epassword", self.__password)

        self.find_and_click(".form_button_submit")  # press login button

        search_button = self.find_elem(".sys-actionbar-action", many=True)[1]
        self.click_elem(search_button)
        # костыль: подождем, пока загрузится через отслеживание
        self.find_elem('.sys-event-hook')
        self.total = self.get_total_number()

    def _paste_dates(self, subm_from: str = None,
                    subm_to: str = None, css=".hasDatepicker"):
        """ .clear() заставляет обновляться таблицу, будем стирать дату
        и мгновенно вставлять новую - так что-то там остается в памяти и все быстро
        """
        fields = self.find_elem(css, many=True)
        if subm_to is not None:
            fields[-1].send_keys(Keys.BACK_SPACE * 10)
            self.find_and_fill(
                css_selector=css,
                text=subm_to,
                many=True,
                field_idx=3,
            )
            self.wait_total_change()
        if subm_from is not None:
            fields[-2].send_keys(Keys.BACK_SPACE * 10)
            self.find_and_fill(
                css_selector=css,
                text=subm_from,
                many=True,
                field_idx=2,
            )
            self.wait_total_change()

    def _activate_base_filtration(self):
        """ too long execution"""
        self.find_and_fill(
            css_selector=".sys-event-hook.sys-fi-mark.yui-ac-input",
            text='Human',
            many=True,
            field_idx=1,
        )
        self.wait_total_change()
        
        for idx in [4, 5]: # 6 is excess
            self.find_elem('.sys-event-hook', many=True)[idx].click()
            self.wait_total_change()
        
        self.paste_dates(subm_to=self.cur_date)
        self.paste_dates(subm_from=self.cur_date)

    def wait_total_change(self, timeout=70):
        wait = WebDriverWait(self.driver, timeout, poll_frequency=1)
        wait.until(lambda x: self.total != self.get_total_number())
        self.total = self.get_total_number()
    
    def _enter_nxt_day(self) -> str:
        nxt_day = next_day_str(self.cur_date)
        self.cur_date = nxt_day
        
        self.paste_dates(subm_to=nxt_day)
        self.paste_dates(subm_from=nxt_day)
        return nxt_day
    
    def button_click(self, button: str, fulltext_mode=False):
        """ button := {reset, fulltext, select, analysis, download} """
        buttons = ["reset", "fulltext", "select", "analysis", "download"]
        assert button in buttons, f"button {button} isn't found"
        b2id = dict(zip(buttons, range(len(buttons))))
        b = self.find_elem(".sys-form-button", many=True)
        b[b2id[button] - int(fulltext_mode)].click()

    def insert_fulltext_search(self, text: str = None, chank_size=1000):
        css = ".sys-event-hook.sys-fi-mark.sys-form-fi-multiline"
        field = self.find_elem(css)
        if text is None:
            field.clear()
        else:
            splitted = []
            for i in range(0, len(text), chank_size):
                chank = text[i:i + chank_size]
                splitted.append(chank)
                field.send_keys(chank)
            assert ''.join(splitted) == text

    def wait_spinners(self, timeout=90, poll_frequency=1):
        sleep(1)
        wait = WebDriverWait(self.driver, timeout, poll_frequency)
        base_spinner = EC.invisibility_of_element_located((By.CLASS_NAME, 'sys_dt_spinner'))
        bad_spinner = EC.invisibility_of_element_located((By.CLASS_NAME, 'sys_timer_inner'))
        wait.until(base_spinner)
        wait.until(bad_spinner)
        
    def all_checkboxes_click(self):
        table_header = self.find_elem(".yui-dt-hd")
        checkbox = table_header.find_element_by_tag_name("input")
        checkbox.click()

    def download_process(self):
        total = self.get_total_number()
        assert total <= 5000, f"Number of choosen records ({total}) > 5000"
        iframe = self.find_elem(".sys-overlay-style")
        self.driver.switch_to.frame(iframe)
        clickables = self.find_elem(".sys-event-hook", many=True)
        clickables[0].click()  # click radio "Input for the Augur pipeline"
        self.wait_spinners()
        clickables = self.find_elem(".sys-event-hook", many=True) # update frame
        clickables[3].click()  # click Download button

        self.driver.switch_to.default_content()
        self.wait_spinners()


def next_day_str(cur_day):
    cur = list(map(int, cur_day.split('-')))
    nxt = datetime.date(*cur) + datetime.timedelta(days=1)
    return str(nxt)


def _collect_totals():
    """ (He)надо глянуть, сколько в каждый день было загрузок """
    # done <= 2020-04-28
    spider = GisaidParser(False)
    spider.login_and_open_table()
    spider.activate_base_filtration()
    
    cur_date = spider.cur_date
    last_date = str(datetime.date.today())  # today
    while cur_date != last_date:
        cur_date = spider.enter_nxt_day()
        cur_total = spider.get_total_number()
        print(f"{cur_date}: {cur_total}")


def num_from_idx(idx):
    return int(re.search('EPI_ISL_(\d+)', idx).groups()[0])


def accertion_idx(first_idx="EPI_ISL_402119", amount=500, exclude=None):
    exclude = exclude or set()
    if not isinstance(exclude, set):
        exclude = set(exclude)
    fnum = num_from_idx(first_idx)
    collection = []
    for x in range(fnum, fnum + amount):
        idx = f"EPI_ISL_{x}"
        if idx not in exclude:
            collection.append(idx)
        if idx == LAST_IDX:
            break
    nxt_idx = f"EPI_ISL_{x + 1}"
    return ",".join(collection), nxt_idx


def idx_iterator(cur_idx='EPI_ISL_573608', amount=500):
    last_idx_num = num_from_idx(LAST_IDX)
    while num_from_idx(cur_idx) < last_idx_num:
        indexes, nxt_idx = accertion_idx(cur_idx, amount)
        yield indexes, cur_idx
        cur_idx = nxt_idx


def start_spider():
    spider = GisaidParser(False)
    spider.login_and_open_table()
    spider.button_click('fulltext')
    sleep(1)
    return spider


def parsing_step(spider, indexes):
    spider.insert_fulltext_search(indexes, INSERT_CHANK_SIZE)
    # spider.wait_total_change()
    spider.wait_spinners()
    sleep(5)
    
    spider.all_checkboxes_click()
    sleep(2)
    spider.wait_spinners()
    
    spider.button_click('download', fulltext_mode=True)
    spider.wait_spinners()
    
    spider.download_process()
    sleep(1)
    spider.all_checkboxes_click()
    spider.wait_spinners()
    spider.insert_fulltext_search()  # clear fulltext field
    spider.wait_spinners()


def main():
    start_idx = read_last_idx()
    spider = start_spider()
    indexes_loader = idx_iterator(start_idx, SEQ_AMOUNT)

    for indexes, cur_idx in indexes_loader:
        write_last_idx(cur_idx)
        while True:
            try:
                parsing_step(spider, indexes)
                print(f"Done: {cur_idx}")
                break
            except Exception as e:
                print(f"Failed: {cur_idx}")
                print(e)
                print("Restart spider\n")
                spider.close()
                del spider
                spider = start_spider()



def meta_parsing_step(spiders, spid, indexes, cur_idx):
    spider = spiders[spid]
    while True:
        try:
            parsing_step(spider, indexes)
            print(f"Done: {cur_idx}")
            break
        except Exception as e:
            print(f"Failed: {cur_idx}")
            print(e)
            print("Restart spider\n")
            spider.close()
            del spider
            spider = start_spider()
            spiders[spid] = spider  #change list of spiders


def main_milti(n_threads: int):
    """http://python-3.ru/page/multiprocessing"""
    start_idx = START_IDX
    indexes_loader = idx_iterator(start_idx, 500)
    spiders = [start_spider() for _ in range(n_threads)]
    while True:
        processes = []
        for spid in range(n_threads):
            indexes, _cur_idx = next(indexes_loader)
            proc = Process(
                target=meta_parsing_step, 
                args=(spiders, spid, indexes, _cur_idx),
            )
            processes.append(proc)
            proc.start()
        
        for proc in processes:
            proc.join()
        

if __name__ == "__main__":
    main()
