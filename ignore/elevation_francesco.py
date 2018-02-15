from pylab import *
import pyrap.quanta as qa
import pyrap.tables as pt
import pyrap.measures as pm
import sys, glob
import numpy

for obs in sorted(glob.glob('c09-o0*/3c380/')):
    msname = glob.glob(obs+'/*MS')[0]
    print msname

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
    if ra<0: ra += 2*numpy.pi
    dec = direction[ ant_no, field_no, 1 ]
    target = {'name' : 'Pointing', 'ra' : ra, 'dec' : dec}
    #print "Target ra/dec (deg):", target['ra']*180/numpy.pi, target['dec']*180/numpy.pi
    field_table.close()

    # Get a ordered list of unique time stamps from the measurement set
    time_table = pt.taql('select TIME from $1 orderby distinct TIME', tables = [ms])
    time = time_table.getcol('TIME')
    time1 = time/3600.0
    time1 = time1 - floor(time1[0]/24)*24

    ra_qa  = qa.quantity( target['ra'], 'rad' )
    dec_qa = qa.quantity( target['dec'], 'rad' )
    pointing =  me.direction('j2000', ra_qa, dec_qa)

    t = qa.quantity(time[0], 's')
    t1 = me.epoch('utc', t)
    me.doframe(t1)

    # Loop through all time stamps and calculate the elevation of the pointing
    time_init = time[1800/4]
    #for t in [time_init,time_init+3600,time_init+2*3600,time_init+3*3600,time_init+4*3600]:
    for t in [time_init,time_init+3600,time_init+2*3600]:
        t_qa = qa.quantity(t, 's')
        t1 = me.epoch('utc', t_qa)
        me.doframe(t1)
        a = me.measure(pointing, 'azel')
        elevation = a['m1']
        print elevation['value']/pi*180
