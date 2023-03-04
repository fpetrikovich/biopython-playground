
def exit_with_error(custom_msg = None, err = None):
    if custom_msg != None:
        print("[ERROR] %s" % (custom_msg))
    if err != None:
        print("[ERROR][MESSAGE] %s" % str(err))
    
    exit(1)


def print_warning(custom_msg):
    print("[Warning] %s\n" % (custom_msg))