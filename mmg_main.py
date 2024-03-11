import turning_test, zig_zag_test

__version__ = "1.0.0"
def main():
    print("Version: ", __version__)
    print(
"**************************\n\
1. turning_test\n\
2. zig_zag_test\n\
**************************")
    test_name = input("Enter the test number: ")
    if test_name == "1":
        turning_test.main()
    elif test_name == "2":
        zig_zag_test.main()
    else:
        print("Invalid test name")

if __name__ == "__main__":
    main()