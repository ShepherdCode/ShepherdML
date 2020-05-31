import sys

infile=sys.argv[1]
with open("part1."+infile,"w") as files1,\
        open("part2."+infile,"w") as files2,\
        open("part3."+infile,"w") as files3,\
        open("part4."+infile,"w") as files4,\
        open("part5."+infile,"w") as files5:

    myfiles=[]
    myfiles.append(sys.stdout)
    myfiles.append(files1)
    myfiles.append(files2)
    myfiles.append(files3)
    myfiles.append(files4)
    myfiles.append(files5)

    next=0
    with open(infile,"r") as f:
        for rawline in f:
            line=rawline.rstrip('\n')
            if line[0]=='>':
                next += 1
                if next > 5:
                    next = 1
                myfile=myfiles[next]
                print(line,file=myfile)
            else:
                print(line,file=myfile)

