import os
import csv
import time
import glob
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.firefox.options import Options as FirefoxOptions

def setup_driver():
    options = FirefoxOptions()
    options.add_argument("--headless")  # Run in background
    download_dir = os.path.abspath("./data/clinvar")
    
    # Ensure download directory exists
    os.makedirs(download_dir, exist_ok=True)
    
    # Configure download preferences
    options.set_preference("browser.download.folderList", 2)
    options.set_preference("browser.download.dir", download_dir)
    options.set_preference("browser.helperApps.neverAsk.saveToDisk", "text/plain")
    
    return webdriver.Firefox(options=options)

def get_clinvar_data(mim_number, driver):
    download_dir = os.path.abspath("./data/clinvar")
    target_file = os.path.join(download_dir, f"{mim_number}.txt")
    
    # Delete existing clinvar_result files to prevent numbering
    for f in glob.glob(os.path.join(download_dir, "clinvar_result*")):
        try:
            os.remove(f)
        except Exception as e:
            print(f"Couldn't delete {f}: {e}")

    try:
        # Navigate to ClinVar
        driver.get(f"https://www.ncbi.nlm.nih.gov/clinvar?term={mim_number}[MIM]")
        
        # Click "Download" -> "Create File"
        WebDriverWait(driver, 20).until(
            EC.element_to_be_clickable((By.XPATH, "//a[text()='Download']"))
        ).click()
        WebDriverWait(driver, 20).until(
            EC.element_to_be_clickable((By.XPATH, "//button[text()='Create File']"))
        ).click()
        
        # Wait for download to complete (check for 60 seconds)
        start_time = time.time()
        max_wait = 60
        downloaded = False
        
        while time.time() - start_time < max_wait:
            # Check for new clinvar_result.txt (ignore .part files)
            clinvar_files = glob.glob(os.path.join(download_dir, "clinvar_result*.txt"))
            clinvar_files = [f for f in clinvar_files if not f.endswith('.part')]
            
            if clinvar_files:
                newest_file = max(clinvar_files, key=os.path.getctime)
                try:
                    os.rename(newest_file, target_file)
                    print(f"Renamed: {os.path.basename(newest_file)} → {mim_number}.txt")
                    downloaded = True
                    break
                except Exception as e:
                    print(f"Retrying rename due to: {e}")
                    time.sleep(2)
            time.sleep(1)
        
        if not downloaded:
            print(f"Timeout for {mim_number}.txt")
            return False
        return True

    except Exception as e:
        print(f"Error for {mim_number}: {e}")
        return False

def main():
    with open('./data/omim_ids.csv', 'r') as f:
        omim_ids = [row["mimNumber"] for row in csv.DictReader(f) if row["mimNumber"] != "Not Found"]
    
    driver = setup_driver()
    
    try:
        for mim in omim_ids:
            output_path = os.path.join("./data/clinvar", f"{mim}.txt")
            if os.path.exists(output_path):
                print(f"Skipping existing: {mim}")
                continue
            
            print(f"Downloading: {mim}")
            success = get_clinvar_data(mim, driver)
            print(f"Status: {'Success' if success else 'Failed'}")
            time.sleep(5)  # Avoid server overload
    
    finally:
        driver.quit()

if __name__ == "__main__":
    main()