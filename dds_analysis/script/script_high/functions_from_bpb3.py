#here are functions copied from bpb3 package
import os

def check_folder(out_folder):
 if not os.path.exists(out_folder):
    print("Create , ", out_folder)
    os.makedirs(out_folder)
 else:
    print("Exists , ", out_folder)

def is_non_zero_file(fpath):
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0


