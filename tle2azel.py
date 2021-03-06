#!/usr/local/bin/python3

'''
Script designed to produce a text output file that will allow <ANTENNA> to track <SATELLITE>

This requires at least one input file in directory called {antenna.input}. Format of 
{antenna.input} is 
!   NAME----  LONG------  LAT------  Z---  AZR  ELR
    HOBART26  +147.44052  -42.80358    65   40   40 

Where columns are [1] 6char antenna name; [2] antenna longtitude (deg) to 5dp inc east direction;
[3] antenna latitude (deg) to 5dp inc direction N|S; [4] height in m above sea level, [5] azimuthal 
max slew rate (deg/min); [6] antenna elevation max slew rate (deg/min)

Other ideas I have is adding:

Plotting - Done 15/02/2021
Single or multiple tracks -n 
Beam size considerations
Transmission frequency of satellites vs. rxp's at antennas
Slew rate considerations vs. satellite speed (e.g. ISS)
Max/min azimuthal limits/wraps
XYZ geodetic coordinates and drift rates

Lucas Jordan Hyland

2020/02/15

'''

###################################################################################################
import datetime, os, sys, argparse, numpy as np
import matplotlib.pyplot as plt, matplotlib.dates as dates
from matplotlib import rc
from astropy import time as astrot, units as u
from astropy.coordinates import EarthLocation
from skyfield.api import Topos, load, EarthSatellite
import pdb

rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
rc('text', usetex=True)

###################################################################################################

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose",
                        help="increase output verbosity",
                        action="count",default=0)
    parser.add_argument('antenna',
                        help='antenna found in antennas.input',
                        type=str)
    parser.add_argument('satellite',
                        help='satellite or spacecraft name or ID number or "up"',
                        type=str)
    parser.add_argument('time',
                        help='start time of observations or "now" \\ yyyy-mm-ddTHH:MM:SS.S',
                        type=str)
    parser.add_argument('-f','--file',
                        help='input file for satellites',
                        type=str, default=None)
    parser.add_argument('-p','--plot',
                        help='plot output flag',
                        action='count',default=0)
    parser.add_argument('-n','--tracks',
                        help='number of tracks from start. n=0 indicates all',
                        action='count',default=0)
    parser.add_argument('-L','--limit',
                        help='antenna elevation limit',
                        type=float,default=10.)
    parser.add_argument('-o','--outfile',
                        help='output file for tracking',
                        type=str,default='tle2azel')
    args = parser.parse_args()
    ###
    global ts
    # this is to prevent ut0-ut1 loading error
    ts = load.timescale(builtin=True) 
    # identify the telescope amd parameters
    name, longi, lati, z, max_az, max_el = _match_antenna(args)
    try:
        ant_ = Topos(latitude=float(lati), longitude=float(longi), elevation_m=float(z))
    except ValueError:
        ant_ = Topos(latitude=lati, longitude=longi, elevation_m=float(z))
    # make a time scale for a 24hour period beginning at given time every 10sec
    t_r = _generate_times(args)
    # first want to load TLEs. If an input file is given, use that. Else look for default TLEs.
    if args.file!=None: satellite_tles = load.tle(args.file)
    else:               satellite_tles = _download_satellite_file()
    # check if satellite id/name given, else print all UP at given time
    # match satellite name and create TLE object
    if args.satellite in ("0","up"): _which_satellite_up(args,ant_,satellite_tles)  
    else: sat_ = _match_satellite(args,satellite_tles)
    # make track with skyfield.api
    track = (sat_ - ant_).at(ts.from_astropy(t_r)).altaz()
    t_jd, t_mjd, az, el, rng = _process_track(args,t_r,track)
    # print out this information into 2 files, HMI format and sattrack format
    _hmi_outprint(args, az, el, rng, t_mjd)
    _sattrakk_outprint(args, az, el, rng, t_jd)
    if args.plot>0: _plottrack(args,t_mjd,az,el)
    print('Processed track for {0:20s} starting {1:30s}'.format(sat_.name,args.time))

###################################################################################################

def _download_satellite_file():
    url     = 'http://celestrak.com/NORAD/elements'
    satfiles= ['gps-ops.txt','gnss.txt','weather.txt','glo-ops.txt']
    infile  = 'satellites.input'
    '''
    Check if satellite input file already exists, if not check if downloadable 
    satellite file exists. If not download and append
    '''
    if os.path.exists(infile):
        os.remove(infile)
    os.system('touch {0:}'.format(infile))
    for sfile in satfiles:
        if not os.path.exists(sfile):
            os.popen('wget --output-document={1:s} {0:s}/{1:s}'.format(url,sfile))
            os.system('cat {0:s} >> {1:s}'.format(sfile,infile))
        else:
            os.system('cat {0:s} >> {1:s}'.format(sfile,infile))
    return load.tle(infile)

def _match_satellite(args,satellite_tles):
    ids, names = [], []
    for item in list(satellite_tles.keys()):
        if type(item)==int:
            ids.append(item)
            names.append(satellite_tles[item].name)
    # if TLE input file only has one satellite, use that and ignore input name.
    # else try and match (should use a while loop here and search all online places)
    if len(names)==1: 
        satellite_name = satellite_names[0]
        satellite      = satellite_tles[satellite_name]
    else:
        try: 
            satellite = satellite_tles[int(args.satellite)]
            return satellite
        except ValueError:
            try:
                satellite = satellite_tles[args.satellite]
                return satellite
            except KeyError:
                print('Invalid satellite identifier {}'.format(args.satellite))
                print('To see more options, increase verbosity (-v)')
                if args.verbose>=1:
                    print('Valid identifiers for your input files are:')
                    print('{0:>6s} : {1:<10s}'.format('ID','NAME'))
                    for j in range(len(ids)):
                        print('{1:>6.0f} : {0:<30s}'.format(names[j],ids[j]))
                sys.exit()

def _match_antenna(args):
    antenna   = 0.
    if not os.path.exists('antennas.input'): 
        sys.exit('Cannot find antennas.input file.')
    for line in open('antennas.input'): 
        if not line.strip().startswith("!"): 
            if line.split()[0]==(args.antenna).upper():
                antenna = line.split()[0]
                return line.split()
    if antenna==0.:
        print('Invalid antenna identifier {}'.format(args.antenna))
        print('To see more options, increase verbosity (-v)')
        if args.verbose>=1:
            print('Valid identifiers are:')
            for line in open('antennas.input'): 
                if not line.strip().startswith("!"): 
                    print('{0:s}'.format(line.split()[0]))
        sys.exit()

def _generate_times(args):
    '''
    Generates observational period or 24hours from start every 10seconds.
    '''
    if not (args.time).lower()=='now': 
        try:
            start = astrot.Time(args.time,format='isot', scale='utc')
        except ValueError:
            print('Invalid time or format {}'.format(args.time))
            print('Valid formats are of the form yyyy-mm-ddTHH:MM:SS.S')
            print('e.g 1989-06-04 or 2012-09-11T04:20:00')
            sys.exit()
    else: start = astrot.Time.now()
    sec_day   = 86400.
    dsec      = 10.
    return start + np.arange(int(sec_day/dsec))*astrot.TimeDelta(dsec, format='sec')

def _hmi_outprint(args,az,el,rng,t):
    if os.path.exists("{0:s}.hmi".format(args.outfile)):
        os.remove("{0:s}.hmi".format(args.outfile))
    #HMI format
    #-862789 091840  0       57513   22604294
    with open("{0:s}.hmi".format(args.outfile), "a") as f:
        print(len(az),file=f)
        for i in range(len(az)):
            print('{0:<+08.0f} {1:07.0f}  {4:20.2f}       {2:5.0f}   {3:8.0f}'.format(az[i]*10000,
                el[i]*10000,int(t[i]),
                24*60*60*1000*( t[i] % int(t[i])),rng[i]),file=f)

def _sattrakk_outprint(args,az,el,rng,t):
    if os.path.exists("{0:s}.sattrak".format(args.outfile)):
        os.remove("{0:s}.sattrak".format(args.outfile))
    #printing sattrak output
    #JD : 2455352.32379623 ; Az : -83.882 ; El : 13.300 ; Rn :  69420.0
    with open("{0:s}.sattrak".format(args.outfile), "a") as f:
        for i in range(len(az)):
            print('JD : {0:16.8f} ; Az :{1:>+8.3f} ; El : {2:6.3f} ; Rn :{3:10.2f}'.format(t[i],
                az[i],el[i],rng[i]),file=f)

def _process_track(args,t_r,track,max_az=40,max_el=40):
    '''
    Convert track object and times into az, el and remove limits and large slew rates
    One day I'll add the slew rate and telescope azimuthal limits here
    '''
    # turn track into time, az, el
    el    = track[0].degrees
    az    = track[1].degrees
    d     = track[2].km
    az_l  = az[el>args.limit]
    el_l  = el[el>args.limit]
    t_mjd = t_r.mjd[el>args.limit]
    t_jd  = t_r.jd[el>args.limit]
    rng   = d[el>args.limit]
    limit = args.limit
    #dazdt = np.diff(az)/np.diff(T)/60. # deg/min
    #deldt = np.diff(al)/np.diff(T)/60. # deg/min
    return t_jd, t_mjd, az_l, el_l, rng
 
def _which_satellite_up(args,antenna,satellites):
    if not (args.time).lower()=='now': 
        try:
            asT = astrot.Time(args.time,format='isot',scale='utc')
        except ValueError:
            print('Invalid time or format {}'.format(args.time))
            print('Valid formats are of the form yyyy-mm-ddTHH:MM:SS.S')
            print('e.g 1989-06-04 or 2012-09-11T04:20:00')
            sys.exit()
    else: asT = astrot.Time.now()
    t  = ts.from_astropy(asT)
    t2 = ts.from_astropy(asT+5*u.min)
    print('{0:<20s} {1:>8s} {2:>10s} {3:>10s}'.format('Name','ID','AZ','R|S'))
    for sat in satellites:
        try: 
            int(sat)
            skip=0
        except ValueError: 
            skip=1
        if skip==1: continue
        topo_t = (satellites[sat] - antenna).at(t).altaz()
        az1 = topo_t[0].degrees
        az2 = (satellites[sat] - antenna).at(t2).altaz()[0].degrees
        if az2>az1: rate='RISING'
        elif az2<az1: rate='SETTING'
        if topo_t[0].degrees<args.limit:continue
        #print(sat)
        print('{0:<30s} {1:>8.0f} {2:>10.2f} {3:>10s}'.format(satellites[sat].name,
            int(sat),topo_t[0].degrees,rate))
    sys.exit()

def _plottrack(args,t,az,el):
    for postfix in ['trak','azel','rate']:
        plotfile = '{}_{}.pdf'.format(args.outfile,postfix)
        if os.path.exists(plotfile): os.remove(plotfile)
    
    fig1 = plt.figure(1,figsize=(5,5))
    ax1  = fig1.add_subplot(111,projection='polar')
    ax1.scatter(az[0]/57. ,el[0] ,marker=',',color='b',s=50) #start
    ax1.scatter(az[-1]/57.,el[-1],marker=',',color='r',s=50) #stop
    ax1.scatter(az/57.    ,el    ,marker='.',color='k',s=1)
    ax1.set_ylim(90,0)
    ax1.set_theta_zero_location('N')
    ax1.set_theta_direction(-1)
    ax1.legend([r'Start',r'Stop'],loc='right')
    fig1.savefig('{}_trak.pdf'.format(args.outfile),bbox_inches='tight')
    ####
    fig2, ax2 = plt.subplots(1,2,figsize=(10,5))
    ax2[0].plot(t,az,'b,')
    ax2[1].plot(t,el,'r,')
    ax2[1].set_ylim([0,90])
    ax2[0].set_ylabel(r'Azimuth (deg)')
    ax2[1].set_ylabel(r'Elevation (deg)')
    fig2.savefig('{}_azel.pdf'.format(args.outfile),bbox_inches='tight')
    ####
    dazdt = np.diff(az)/np.diff(t)/(24.*60) # deg/min
    deldt = np.diff(el)/np.diff(t)/(24.*60) # deg/min
    t_d   = t[:-1] + np.diff(t) 
    #
    fig3, ax3 = plt.subplots(1,2,figsize=(10,5))
    ax3[0].plot(t_d,dazdt,'b,')
    ax3[1].plot(t_d,deldt,'r,')
    ax3[0].set_ylabel(r'Az Rate (deg/min)')
    ax3[1].set_ylabel(r'El Rate (deg/min)')
    for AX in ax3:
        AX.set_ylim([-40,40])
    fig3.savefig('{}_rate.pdf'.format(args.outfile),bbox_inches='tight')


###################################################################################################

if __name__=='__main__':
    main()

###################################################################################################