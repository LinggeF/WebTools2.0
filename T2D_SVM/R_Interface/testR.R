
args=commandArgs(TRUE);

n_modules=args[1]
ga_times=args[2]
ls_times=args[3]
input1=args[9]
input2=args[10]
input3=args[11]
output_path=args[12]

print("n_modules:")
print(n_modules)
print("ga_times:")
print(ga_times)
print("ls_times:")
print(ls_times)
print("input1 path:")
print(input1)
print("input2 path:")
print(input2)
print("input3 path:")
print(input3)
print("output_path")
print(output_path)

lines1=readLines(input1)
lines2=readLines(input2)
lines3=readLines(input3)
os=file(output_path,"w")
write(n_modules,os,append = T)
write(ga_times,os,append = T)
write(ls_times,os,append = T)
write(lines1,os,append = T)
write(lines2,os,append = T)
write(lines3,os,append = T)
close(os)
