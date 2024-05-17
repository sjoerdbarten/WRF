diri = "./"
flnm = "Vredepeel-NOx.csv"
strs = asciiread(diri+flnm,-1,"string")
delim = ","
concs = str_get_field(strs,2,delim)

print(concs)
