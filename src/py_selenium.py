
codes = '''
import time

from selenium import webdriver
from selenium.webdriver.common.action_chains import ActionChains
import pyautogui

driver = webdriver.Chrome()
driver.get("https://data.sdss.org/sas/dr12/boss/")
ActionChains(driver).click(driver.find_element_by_link_text('lss/')).perform()
ActionChains(driver).context_click(driver.find_element_by_link_text('boss_lss.sha1sum')).perform()
#time.sleep(5)
pyautogui.typewrite(['down', 'down', 'down', 'down'])
time.sleep(1)
pyautogui.typewrite(['return'])
#time.sleep(1)
pyautogui.moveTo(x=1400, y=965,duration=3, tween=pyautogui.linear)
pyautogui.click()
#time.sleep(5)
#ActionChains(driver).context_click(driver.find_element_by_link_text('Save')).perform()
'''

print(codes)
