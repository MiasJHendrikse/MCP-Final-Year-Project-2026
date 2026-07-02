# Open with UTF-8 and ignore decoding hiccups
with open("uae6.z18.00.s0500000.eng", "r", encoding="utf-8", errors="ignore") as file:
    for i in range(20):
        print(file.readline())