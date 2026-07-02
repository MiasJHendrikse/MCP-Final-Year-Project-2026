import os

_data_path = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "data", "uae6.z18.00.s0500000.eng",
)

# Open with UTF-8 and ignore decoding hiccups
with open(_data_path, "r", encoding="utf-8", errors="ignore") as file:
    for i in range(20):
        print(file.readline().encode("ascii", errors="replace").decode("ascii"))
