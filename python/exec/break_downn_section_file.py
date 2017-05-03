#!/usr/bin/env python

#       B a r a K u d a
#
#
#
#       L. Brodeau, 2017
#
import sys
from os import path

csn = sys.argv[0]

def get_sections(cfile):
    list_sections = []
    list_coord    = []
    list_refts    = []
    f = open(cfile, 'r') ; cread_lines = f.readlines() ; f.close()
    jl=0 ; l_stop=False
    while not l_stop:        
        ll = cread_lines[jl]; ls = ll.split(); cc = ls[0]
        if cc == 'EOF':
            l_stop=True            
        elif cc[0] != '#':
            if cc != 'EOF' and cc[0:13] != 'ref_temp_sali':
                #print cc
                list_sections.append(cc)
                #
                # now time to read the coordinates:
                jl=jl+1
                lok=False
                while not lok:
                    ll = cread_lines[jl]; ls = ll.split(); cc = ls[0]
                    #print ' ls => ', ls
                    #print ' cc => ', cc
                    if cc == 'EOF' or cc[0:13] == 'ref_temp_sali':                 
                        print 'break_downn_section_file.py get_sections ERROR1!'; sys.exit(0)
                    elif cc == '#':
                        print '#'
                    else:
                        #print ' to append ls => ', ls
                        # Read
                        list_coord.append(ls)
                        lok=True
                        # Check if ref_temp_sali to read:
                        ll = cread_lines[jl+1]; ls = ll.split(); cc = ls[0]
                        if cc[0:13] == 'ref_temp_sali':
                            list_refts.append(ls)
                            jl=jl+1
                        else:
                            list_refts.append(['no','no','no'])

        else:
            print '  ....  '


        jl=jl+1    
    return list_sections, list_coord, list_refts





narg = len(sys.argv)
if narg != 3:
    print 'Usage: {} <ASCII section file> <max_number_sect_per_file>'.format(csn)
    sys.exit(0)
cf_in = sys.argv[1]
cdum  = sys.argv[2]

max_number_sect_per_file = int(cdum)


if not path.exists(cf_in): print ' File: '+cf_in+' is not there!'; sys.exit(0)


list_sn, list_co, list_ts = get_sections(cf_in)

#list_sections = bt.get_sections_from_file(cf_in)


print ''

#print '\n\n Section:\n', list_sn[:]
#print '\n\n Coor:\n', list_co[:]
#print '\n\n Ref:\n', list_ts[:]

nbs = len(list_sn)
#print ' *** Number of sections                  =', nbs
#print ' *** Max number of sections in new files =', max_number_sect_per_file

nbfiles = nbs/max_number_sect_per_file
if nbfiles*max_number_sect_per_file < nbs: nbfiles = nbfiles+1
#print ' *** nbfiles =', nbfiles ; # number of sections per file!



#lok = False
#while not lok:
#nbfiles = nbs/(max_number_sect_per_file)
#if nbs%(max_number_sect_per_file) > nbfiles: nbfiles = nbfiles+1



#n_extra_last_file = 0
#if nbfiles*max_number_sect_per_file < nbs: n_extra_last_file = nbs - nbfiles*max_number_sect_per_file
#print ' *** n_extra_last_file =', n_extra_last_file



jcum = 0

for jf in range(nbfiles):
    clab = '%2.2i'%(jf+1)

    cf_out = 'transportiz_'+clab+'.dat'
    f = open(cf_out, 'wb')

    for jsf in range(max_number_sect_per_file):
        if jcum < nbs:
            f.write(list_sn[jcum]) ; f.write('\n')

            vv = list_co[jcum] ; cs=''
            for cv in vv:
                f.write(cs+cv); cs=' '
            f.write('\n')    

            vv = list_ts[jcum] ; cs=''
            if vv[0] != 'no':
                for cv in vv: f.write(cs+cv); cs=' '
                f.write('\n')
        jcum = jcum + 1
    
    f.write('EOF\n')
    f.close()
    print ' *** '+cf_out+' written!\n'






sys.exit(0)

f = open('data_new.txt', 'wb')
for js in range(nbs):
    f.write(list_sn[js]) ; f.write('\n')
    vv = list_co[js] ; cs=''
    for cv in vv:
        f.write(cs+cv); cs=' '
    f.write('\n')

    vv = list_ts[js] ; cs=''
    if vv[0] != 'no':
        for cv in vv:
            f.write(cs+cv); cs=' '
        f.write('\n')
            

f.write('EOF\n')
f.close()





















print ''+csn+' done...\n'







#  LocalWords:  jl
