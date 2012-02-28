"""Helper classes for Python functional programming."""

def coroutine(func):
    """A decorator function that eliminates a need to call .next() to start coroutine.
    Example of use:
    @coroutine
    def g():
        y = 0
        x = (yield)
        while 1:
            y += x
            x = (yield y)
    ig = g()
    #calling next() is not needed anymore
    #print ig.next() #prints None
    print ig.send(1) #prints 1
    print ig.send(2) #prints 3
    """
    def start(*args,**kwargs):
        cr = func(*args,**kwargs)
        cr.next()
        return cr
    return start

