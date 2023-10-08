import pandas as pd
import requests
import re
import time
from bs4 import BeautifulSoup
import os
from selenium import webdriver
import time
from selenium.webdriver.common.by import By
import sys

def data_scraping(output_path):
    """Scrape the data from homegate.ch and save it to a csv file.
    
    Args:
        output_path (str): Path to the raw data (Homegate_data.csv)
    """

    chromedriver_path = os.getcwd() + '/chromedriver'
    browser = webdriver.Chrome(executable_path=chromedriver_path)

    # Go through all 1-50 of the search pages and find all the links in each page:
    urls = []
    for i in range(1, 50):    
        search_url_pattern = 'https://www.homegate.ch/rent/real-estate/city-zurich/matching-list?ep={}&be=2000&o=dateCreated-desc'
        search_url = search_url_pattern.format(i)
        response = browser.get(search_url)
        elems = browser.find_elements(By.TAG_NAME, "a")
        for elem in elems:
            urls.append(elem.get_attribute("href"))
    
    # Narrow the links down to only those that are of rental property listings:
    pattern = 'https://www.homegate.ch/rent/\d'
    urls2 = [s for s in urls if re.match(pattern, s)]
    df = pd.DataFrame(data={"urls": urls2})
    df.to_csv("data/urls_list.csv", sep=',',index=False)

    # We have 980 urls which we should check with selenium. In the following we scrape the multiple sites. This took about 1.5 hours.
    df2 = pd.DataFrame()
    counter = 0
    for hg_url in urls2:
        try:
            response = browser.get(hg_url)
            time.sleep(5)
            data = [dd.text for dd in browser.find_elements(By.TAG_NAME, 'dd')]
            categories = [dt.text for dt in browser.find_elements(By.TAG_NAME, 'dt')]
            df = pd.DataFrame([dict(zip(categories, data))])
            df['Address'] = browser.find_element(By.CLASS_NAME, 'AddressDetails_addressDetails_wuB1A .AddressDetails_address_i3koO').text
            df['URL'] = hg_url
            df2 = pd.concat([df2, df], ignore_index = True)   
        except: 
            pass
        counter += 1
        if counter%10 == 0:
            print(counter)
            df2.to_csv(output_path, index = False)


if __name__ == "__main__":
    data_scraping(
        snakemake.output[0]
    )