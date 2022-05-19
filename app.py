from re import I
from flask import Flask


app = Flask(__name__)


@app.route("/data")
def data():
    with open("data.csv", "r") as f:
        data = f.read()
    return data


if __name__ == "__main__":
    app.run()
