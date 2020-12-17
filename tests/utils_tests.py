def assert_files_are_equal(file1, file2):
    with open(file1, "r") as f1:
        with open(file2, "r") as f2:
            content_1 = f1.read()
            content_2 = f2.read()
            assert content_1 == content_2