from noaatools.imageproc import process_img
import unittest
import os
import subprocess


class Georefests(unittest.TestCase):
    def setUp(self):
        self.get_test_data()

    TEST_URL = "https://svarog.space/data/ab802ca1-419f-418a-aeea-d99bb9c702aa-0-1c60fdf5-3f18-409d-9834-17f014c608c1_product.png"
    TEST_FILE = "tests/1276.png"

    def get_test_data(self):
        """This ensures the data file is here. We don't want to commit 3MB data to this repo."""
        try:
            os.stat(self.TEST_FILE)
            print(f"File {self.TEST_FILE} exists, skipping download.")
        except FileNotFoundError:
            print("Test data file %s missing, downloading..." % self.TEST_FILE)
            args = ["wget -nd %s -O %s" % (self.TEST_URL, self.TEST_FILE)]
            process = subprocess.run(args, capture_output=True, shell=True)
            if process.returncode != 0:
                print("wget failed to download file. details:\n===STDOUT===")
                print(process.stdout)
                print("===STDERR===")
                print(process.stderr)
            self.assertEqual(process.returncode, 0)
            print("Output file stored in %s" % self.TEST_FILE)

    def test_ellipsoids(self):

        params = {
            "histogram": True,
            "histogram-adaptive": True,
            "border": True,
            "show": False,
            "write": True,
            "write-left": True,
            "write-right": False,
            "denoise": False,
            "georef": True  # Georeference
        }

        process_img("tests/1276.png", params)

        # self.assertIs()
