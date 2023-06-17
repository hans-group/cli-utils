def prog_bar(iterator, total):
    def make_progbar(num_done, total_length):
        num_done_syms = int((num_done + 1) / total_length * 30)
        if total_length > 1 and num_done == 0:
            progbar = "[" + ">".ljust(30, ".") + "]"
        elif num_done < total_length - 1:
            progbar = (
                "[" + "=" * num_done_syms + ">".ljust(30 - num_done_syms, ".") + "]"
            )
        else:
            progbar = "[" + ">".rjust(30, "=") + "]"

        return progbar

    for i, item in enumerate(iterator):
        yield item
        print(make_progbar(i, total), end="\r")
    print(make_progbar(i + 1, total))
