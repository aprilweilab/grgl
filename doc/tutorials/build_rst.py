import glob
import subprocess

def convert_to_rst(notebook: str, rst_directory: str):
    subprocess.check_call(["jupyter", "nbconvert", "--to", "rst", "--output-dir", rst_directory, notebook])

if __name__ == "__main__":
    for notebook in glob.glob("notebooks/*.ipynb"):
        print(f"Converting {notebook} to ReStructured Text")
        convert_to_rst(notebook, "./")
