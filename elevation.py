from pylab import *
#import pyrap.quanta as qa
#import pyrap.tables as pt
#import pyrap.measures as pm
import casacore.quanta as qa
import casacore.tables as pt
import casacore.measures as pm
import sys, glob
import numpy as np

#for obs in sorted(glob.glob('c09-o0*/3c380/')):
msname = sys.argv[1]
print 'MS file: ', msname

# Create a measures object
me = pm.measures()

# Open the measurement set and the antenna and pointing table
ms = pt.table(msname, ack=False)

# Get the position of the first antenna and set it as reference frame
ant_table = pt.table(msname + '/ANTENNA', ack=False)
ant_no = 0
pos = ant_table.getcol('POSITION')
x = qa.quantity( pos[ant_no,0], 'm' )
y = qa.quantity( pos[ant_no,1], 'm' )
z = qa.quantity( pos[ant_no,2], 'm' )
position =  me.position( 'wgs84', x, y, z )
me.doframe( position )
#print position
ant_table.close()

# Get the first pointing of the first antenna
field_table = pt.table(msname + '/FIELD', ack=False)
field_no = 0
direction = field_table.getcol('PHASE_DIR')
ra = direction[ ant_no, field_no, 0 ]
if ra<0: ra += 2*np.pi
dec = direction[ ant_no, field_no, 1 ]
target = {'name' : 'Pointing', 'ra' : ra, 'dec' : dec}
print "Target ra/dec (deg):", target['ra']*180/np.pi, target['dec']*180/np.pi
field_table.close()

# Get a ordered list of unique time stamps from the measurement set.
time_table = pt.taql('select TIME from $1 orderby distinct TIME', tables = [ms])
time = time_table.getcol('TIME')


ra_qa  = qa.quantity( target['ra'], 'rad' )
dec_qa = qa.quantity( target['dec'], 'rad' )
pointing =  me.direction('j2000', ra_qa, dec_qa)

t = qa.quantity(time[0], 's')
t1 = me.epoch('utc', t)
print 't1=', t1
me.doframe(t1)

# Loop through all time stamps and calculate the elevation of the pointing
time_init = time[0]
elev = []
times = np.arange(time[0], time[-1], 100)
for t in times:
    t_qa = qa.quantity(t, 's')
    t1 = me.epoch('utc', t_qa)
    me.doframe(t1)
    a = me.measure(pointing, 'azel')
    elevation = a['m1']
    print elevation['value']/pi*180
    elev.append(elevation['value']/pi*180)

from matplotlib.pyplot import figure, show

# Convert from seconds to hours.
times_hours = times/3600.0
# Convert from hours to days, round down to the closest day, convert back to hours and subtract.
# This gives time of day in hours.
times_hours = times_hours - floor(times_hours[0]/24)*24

fig = figure()
fig.suptitle(msname, fontweight='bold')
ax = fig.add_subplot(111)
ax.scatter(times, elev, marker='.', color='k')
ax.set_xlabel('Time [s]')
ax.set_ylabel('Elevation [deg]')
fig.savefig('elevation.png', bbox_inches='tight', dpi=250)
show()
