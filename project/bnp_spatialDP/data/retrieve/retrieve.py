# http://stackoverflow.com/questions/13260476/how-to-pass-input-to-a-web-page-using-a-automated-script
site="https://eosweb.larc.nasa.gov/cgi-bin/sse/daily.cgi?email=&step=1&lat=666&lon=666&sitelev=&ms=1&ds=1&ys=1985&me=12&de=31&ye=2004&p=swv_dwn&submit=Submit&plot=swv_dwn"
import urllib
import urllib2
import string
import re # Regexpr
import numpy as np
import sys

gridlocs = np.genfromtxt("gridlocs.dat",delimiter=',',skip_header=1)
n = len(gridlocs)

for i in range(0,n):
    lat = gridlocs[i,1]
    lon = gridlocs[i,0]
    values = {'lat':lat, 'lon':lon, 'ms':'1', 'ds':'1', 'ys':'1985', 'me':'12', 'de':'31', 'ye':'2004', 'submit':'Submit', 'param':'TSKIN'}
    data = urllib.urlencode(values)
    req = urllib2.Request(site,data)
    #
    opener1 = urllib2.build_opener()
    page1=opener1.open(req)
    #
    htmlfile=page1.read()
    txt = re.search("<a href=.*Download a text file",htmlfile)
    txtpath = "https://eosweb.larc.nasa.gov" + txt.group(0).split('"')[1]
    txtfile = urllib2.urlopen(txtpath)
    out = txtfile.read()
    #out2 = out.split('\n')
    #out2 = out2[8:len(out2)-1]
    #fout.write('\n'.join(out2))
    fout = open('temps/'+`i`+'.dat', "wb")
    fout.write(out)
    fout.close()
    #print "\rProgress: %d%s" %(round(i*100/n),"%"),
    print "\rProgress: %d%s%d" %(i,"/",n),
    sys.stdout.flush()
