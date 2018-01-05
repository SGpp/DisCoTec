from neovim import attach

import os, tempfile, random, time, sys

def pipe(input, output):
    while 1:
        line = input.readline()
        output.write(line)
        output.flush()

def streamStdInToNvim():
    nvim = attach('socket', path='/tmp/nvim')

    tmpdir = tempfile.mkdtemp()
    filename = os.path.join(tmpdir, 'myfifo' + str(random.randint(0, 2**16)))
    print(filename)
    try:
        os.mkfifo(filename)
    except OSError as e:
        print("Failed to create FIFO: %s" % e)
        return
    try:
        # print('vsplit ' + filename)
        cmd = r'tabedit term://bash --rcfile <(echo \"cat ' + filename + r'\" )'
        nvim.command(cmd)
        with open(filename, 'w') as fifo:
            pipe(sys.stdin, fifo)
    finally:
        os.remove(filename)
        os.rmdir(tmpdir)

if __name__ == "__main__":
    streamStdInToNvim()
