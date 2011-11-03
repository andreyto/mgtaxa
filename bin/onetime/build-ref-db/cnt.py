import fileinput
cnt = set()
for l in fileinput.input():
    cnt.add(l.split()[1])
print sorted(cnt)

