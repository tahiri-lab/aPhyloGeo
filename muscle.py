import subprocess

def fetchingData():
    subprocess.call("./fetch_data.sh")

if __name__ == '__main__':
    fetchingData()