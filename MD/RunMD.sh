#!/bin/bash
# rev_prefix=`svn info |grep Revision |awk 'print $2'`
svn update
rev_prefix=`svnversion`
exec=MD
# Usage: MD [prefix] [chain] [temperature] [total_time] [dt] [epsi] [q] [Ec] [print_each]
chain=1000100010001
temperature=0.04
total_time=20000
dt=0.001
epsi=1.0
q=0.1
Ec=-1.0
print_each=100
${exec} ${rev_prefix} ${chain} ${temperature} ${total_time} ${dt} ${epsi} ${q} ${Ec} ${print_each}


