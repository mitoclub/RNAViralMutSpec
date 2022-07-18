from time import sleep

from selenium.webdriver.firefox.options import Options as FirefoxOptions
from selenium.webdriver.chrome.options import Options as ChromeOptions
from selenium.webdriver.support.ui import WebDriverWait
from selenium import webdriver

FIREFOXDRIVER_PATH = "./scripts/parser/geckodriver"
CHROMEDRIVER_PATH = "./scripts/parser/chromedriver"


class Parser():
    """
    https://stackoverflow.com/questions/6682503/click-a-button-in-scrapy
    https://stackoverflow.com/a/60627463/14998254
    https://selenium-python.readthedocs.io/waits.html
    """

    def __init__(self, headless, wdriver='chrome', driver_path=CHROMEDRIVER_PATH):
        self.driver = self.run_driver(
            headless=headless,
            executable_path=driver_path,
            wd=wdriver,
        )

    @staticmethod
    def run_driver(headless: bool, executable_path: str, wd: str):
        if wd == 'firefox':
            options = FirefoxOptions()
            dclass = webdriver.Firefox
        elif wd == 'chrome':
            options = ChromeOptions()
            dclass = webdriver.Chrome
        else:
            raise ValueError(f"Incorrect driver type '{wd}'")
        
        if headless:
            options.add_argument("--headless")
        
        driver = dclass(options=options, executable_path=executable_path)
        return driver

    @staticmethod
    def _run_driver(headless: bool, executable_path: str, wd: str):
        if wd == 'firefox':
            options = FirefoxOptions()
            if headless:
                options.add_argument("--headless")
            dclass = webdriver.Firefox
            return dclass(options=options,executable_path=executable_path)

        elif wd == 'chrome':
            return webdriver.Chrome(executable_path=executable_path)
        else:
            raise Exception(f'There are no webdriver {wd}')
            
        driver = dclass(options=options,executable_path=executable_path)
        return driver

    def click_elem(self, elem, time=0):
        if time == 5:
            raise TimeoutError(
                f'Cannot click button-element after {time} attempts')
        try:
            elem.click()
        except BaseException:
            sleep(5)
            self.click_elem(elem, time + 1)

    def _find_elem(self, css_selector, many=False, response=None, time=0):
        """
        params:
            - css_selector
            - many: one or many elements to return
            - response: what object use for searching, default self.driver page
            - time: attempt number of searching (against dynamic page loading)
        """
        if time == 5:
            raise TimeoutError(
                f'Cannot find element by css selector after {time} attempts')
        try:
            response = response or self.driver
            if many:
                elem = response.find_elements_by_css_selector(css_selector)
            else:
                elem = response.find_element_by_css_selector(css_selector)
            return elem
        except BaseException:
            sleep(5)
            return self._find_elem(css_selector, many, response, time + 1)

    def find_elem(self, css_selector, many=False, response=None):
        response = response or self.driver
        wait = WebDriverWait(response, 10)
        if many:
            elem = wait.until(lambda x: x.find_elements_by_css_selector(css_selector))
        else:
            elem = wait.until(lambda x: x.find_element_by_css_selector(css_selector))
        return elem

    def find_and_click(self, css_selector: str):
        self.click_elem(self.find_elem(css_selector))

    def find_and_fill(self, css_selector: str, text: str,
                      many=False, field_idx: int = None):
        """"""
        field = self.find_elem(css_selector, many=many)
        if many:
            field = field[field_idx]
        field.send_keys(text)
    
    def close(self):
        self.driver.close()
