def windowSizePrompt():
    while True:
        size = input("Sliding window size:\n")
        try:
            size = int(size)
        except Exception:
            print("The window size must be a number.\n")
            continue
        if(size <= 0):
            print("The window size must be between greater than 0.\n")
        else:
            break
    return size

def stepSizePrompt(size):
    while True:
        count = input("Step count:\n")
        try:
            count = int(count)
        except Exception:
            print("The step count must be a number.\n")
            continue
        if(count <= 0):
            print("The window size must be between greater than 0.\n")
        elif(count > int(size)):
            print("The step count should be less than the window size.\n")
        else:
            break
    print("SUMMARY")
    print("The window size is: "+ str(size) +".")
    print("The step count is: " + str(count) + ".")




if __name__ == "__main__":
    window_size = windowSizePrompt()
    stepSizePrompt(window_size)
