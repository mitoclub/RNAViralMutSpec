#!/usr/bin/env python3

import datetime
import os
import re
from time import sleep
from typing import Set, Iterable

import click
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.ui import WebDriverWait

from base import Parser

GISAID_URL = "https://www.gisaid.org/"
ACCESS_FILE_PATH = "./scripts/parser/access_file.txt"
LAST_IDX_PATH = "./scripts/parser/last.id"
PATH_TO_ACCESSION_NUMBERS = "./data/Accession_Numbers_omicron_2022_01_13.csv"

START_IDX = "EPI_ISL_1044108"
# LAST_IDX =  "EPI_ISL_5144811"

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

    def wait_spinners(self, timeout=120, poll_frequency=1):
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

    def download_process(self, augur=True):
        total = self.get_total_number()
        max_amount = 5000 if augur else 10000
        assert total <= max_amount, f"Number of choosen records ({total}) > {max_amount}"
        iframe = self.find_elem(".sys-overlay-style")
        self.driver.switch_to.frame(iframe)
        clickables = self.find_elem(".sys-event-hook", many=True)
        if augur:
            clickables[0].click()  # click radio "Input for the Augur pipeline"
            self.wait_spinners()
            clickables = self.find_elem(".sys-event-hook", many=True) # update frame
            clickables[3].click()  # click Download button
        else:
            clickables[4].click()  # click Download button

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


def form_accertion_indexes(first_idx: str, stop_idx: str, amount: int, exclude: Iterable = None):
    exclude = exclude or set()
    if not isinstance(exclude, set):
        exclude = set(exclude)

    fnum = num_from_idx(first_idx)
    collection = []
    for x in range(fnum, fnum + amount):
        idx = f"EPI_ISL_{x}"
        if idx not in exclude:
            collection.append(idx)
        if idx == stop_idx:
            break
    nxt_idx = f"EPI_ISL_{x + 1}"
    return ",".join(collection), nxt_idx


def idx_iterator(start_idx, stop_idx, amount: int):
    last_idx_num = num_from_idx(stop_idx)
    while num_from_idx(start_idx) < last_idx_num:
        indexes, nxt_idx = form_accertion_indexes(start_idx, stop_idx, amount)
        yield indexes, start_idx
        start_idx = nxt_idx


def read_accession_numbers(path: str) -> Set[str]:
    indexes = set()
    with open(path) as fin:
        for line in fin:
            indexes.add(line.rstrip())
    return indexes


def idx_iterator_from_file(records: Set[str], amount: int):
    last_idx = read_last_idx()
    pass_first_recs = last_idx in records

    records: list = sorted(records)
    init_rec_id = records.index(last_idx) if pass_first_recs else 0
    start_idx = records[init_rec_id]

    indexes = []
    for i in range(init_rec_id, len(records)):
        indexes.append(records[i])
        if len(indexes) == amount:
            indexes_str = ",".join(indexes)
            yield indexes_str, start_idx
            
            indexes = []
            start_idx = records[i + 1]
    
    indexes_str = ",".join(indexes)
    yield indexes_str, start_idx


def start_spider():
    spider = GisaidParser(False)
    spider.login_and_open_table()
    spider.button_click('fulltext')
    sleep(1)
    return spider


def parsing_step(spider: GisaidParser, indexes, insert_chank_size, augur=True):
    spider.insert_fulltext_search(indexes, insert_chank_size)
    # spider.wait_total_change()
    spider.wait_spinners()
    sleep(5)
    
    spider.all_checkboxes_click()
    sleep(2)
    spider.wait_spinners()
    
    spider.button_click('download', fulltext_mode=True)
    spider.wait_spinners()
    
    spider.download_process(augur)
    sleep(1)
    spider.all_checkboxes_click()
    spider.wait_spinners()
    spider.insert_fulltext_search()  # clear fulltext field
    spider.wait_spinners()
    sleep(11)


@click.command("gisaid-parser")
@click.option("-a", "--accession_numbers", default=None, type=click.Path(exists=True), help="Path to csv file with Accession Numbers")
@click.option("-1", "--start-idx", default=None, help="If need to download sequences without filtration, pass start idx")
@click.option("-2", "--stop-idx", default=None, help="If need to download sequences without filtration, pass stop idx")
@click.option("--seq-amount", default=SEQ_AMOUNT, type=int, show_default=True, help="Number of seqs to download")
@click.option("--insert-chank-size", default=INSERT_CHANK_SIZE, show_default=True, help="Number of characters to insert into fulltext field at once")
@click.option("--augur", is_flag=True, help="Format of output for Augur pipeline (included metadata) [bool]")
def main(
        accession_numbers: str, 
        start_idx: str, stop_idx: str, 
        seq_amount: int, insert_chank_size: int, augur: bool):

    if accession_numbers is not None:
        indexes_collection = read_accession_numbers(accession_numbers)
        indexes_loader = idx_iterator_from_file(indexes_collection, seq_amount)
    else:
        start_idx = start_idx or read_last_idx()
        assert stop_idx is not None, "start and stop index must be in arguments"
        indexes_loader = idx_iterator(start_idx, stop_idx, seq_amount)

    start_new_spider = True
    for indexes, cur_idx in indexes_loader:
        write_last_idx(cur_idx)
        while True:
            try:
                if start_new_spider:
                    spider = start_spider()
                    start_new_spider = False
                parsing_step(spider, indexes, insert_chank_size, augur)
                print(f"Done: {cur_idx}")
                break
            except Exception as e:
                print(f"Failed: {cur_idx}")
                print(e)
                print("Restart spider\n")
                spider.close()
                del spider
                start_new_spider = True


if __name__ == "__main__":
    main()
