import glob
import subprocess
import sys

def convert_to_rst(notebook: str, rst_directory: str):
    subprocess.check_call(["jupyter", "nbconvert", "--to", "rst", "--output-dir", rst_directory, notebook])

if __name__ == "__main__":
    prefix = sys.argv[1] if len(sys.argv) > 1 else ""
    for notebook in glob.glob("notebooks/*.ipynb"):
        if notebook.startswith(prefix):
            print(f"Converting {notebook} to ReStructured Text")
            convert_to_rst(notebook, "./")
