import sys, time
def progress(iterable,steps=None):
    start = time.time()
    if steps is None:
        steps = len(iterable)
    width = 70
    prog = 0.0
    rest = len("0% (0 of {0}) || Time: ".format(steps))
    block = int(round((width-rest)*prog/steps))
    text = "\r0% (0 of {0}) |{1}| Time: 0:00".format(steps,"-"*(width-rest))
    sys.stdout.write(text)
    sys.stdout.flush()
    for event in iterable:
        yield(event)
        prog += 1
        t = int(time.time()-start)
        sec = str(t%60)
        rest = len("{0}% ({1} of {2}) || Time: ".format(int(prog/steps*100),int(prog),steps))
        block = int(round((width-rest)*prog/steps))
        text = "\r{0}% ({1} of {2}) |{3}| Time: {4}:{5}".format(int(prog/steps*100),int(prog),steps,"#"*block + "-"*(width-block-rest),
                                                                t/60,"0"*(2-len(sec))+sec)
        sys.stdout.write(text)
        sys.stdout.flush()
    sys.stdout.write("\r\n")
    sys.stdout.flush()

if __name__ == "__main__":
    for i in progress(range(30)):
        time.sleep(0.6)
    input()
