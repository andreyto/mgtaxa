from MGT.Util import *

def test_path_hasher():
    names = ["a.txt","a/b.txt","a/b/c.txt","d","e","k.txt","l"]
    top_names = ["a.txt","a","d","e","k.txt","l"] 
    glob_names = ["a.txt","k.txt"]
    root = os.path.join(os.getcwd(),"path_hasher")
    ph = PathHasher(root=root,n_subdirs=2,mode="w")
    paths = {}
    print "Hashing paths..."
    for name in names:
        path = ph(name)
        print path
        paths[name] = path
        if not os.path.dirname(name):
            assert os.path.isdir(os.path.dirname(path))
        dpath = os.path.dirname(path)
        if not os.path.exists(dpath):
            os.makedirs(dpath) #create subdirds in 'name'
        with open(path,"w") as out:
            out.write("ok\n")
    ph = PathHasher(root=root,mode="r")
    for name in names:
        path = ph(name)
        assert paths[name] == path
    print
    print "Listdir..."
    for path in ph.listdir():
        print path
    paths_iter = set((p for p in ph.listdir()))
    paths_names = set((ph(n) for n in top_names))
    assert paths_iter == paths_names
    print
    print "Glob..."
    glob_names_res = []
    for path in ph.glob("*.txt"):
        print path
        glob_names_res.append(os.path.basename(path))
    assert glob_names == glob_names_res

if __name__ == "__main__":
    test_path_hasher()

