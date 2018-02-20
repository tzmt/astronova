from ctypes import *
from ctypes.wintypes import *
import math
import os

# https://docs.python.org/3/library/tkinter.html#tkinter-modules
# http://www.tkdocs.com/tutorial/concepts.html
from tkinter import *
from tkinter import ttk
from tkinter import messagebox

# https://github.com/MrMinimal64/timezonefinder
from timezonefinder import TimezoneFinder

# https://pendulum.eustace.io/docs/#installation
import pendulum

# depencies for pendulum
from datetime import datetime, timedelta

# https://pypi.python.org/pypi/geopy
from geopy.geocoders import Nominatim


# type and constant definitions for interaction with the DLL, and HTML frontend
SIDEREALMODE = c_int32(64*1024)
PLANETLIST = ["Sun", "Moon", "Mercury", "Venus", "Mars", 
              "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"]
CAMPANUS = c_int(67)
VERSION_NUMBER = "0.02a"

# class to hold all natal chart data
class Natal:
    def __init__(self, name):
        self.name = name
        self.datetimelocation = {
            "Year" : 0,
            "Month" : 0,
            "Day" : 0,
            "Hour" : 0,
            "Min" : 0,
            "Sec" : 0,
            "Longitude" : (),
            "Latitude" : ()
            }
        self.planet_dictionary = {
            # The order of the doubles is:
            # Ecliptical longitude, celestial latitude, distance, 
            # Speed in long, speed in lat, speed in dist
            "Sun" : (c_double * 6)(),
            "Moon" : (c_double * 6)(),
            "Mercury" : (c_double * 6)(),
            "Venus" : (c_double * 6)(),
            "Mars" : (c_double * 6)(),
            "Jupiter" : (c_double * 6)(),
            "Saturn" : (c_double * 6)(),
            "Uranus" : (c_double * 6)(),
            "Neptune" : (c_double * 6)(),
            "Pluto" : (c_double * 6)(),
            "Obliquity" : (c_double * 6)(),
            "SVP" : (),
            "LST" : (),
        }
        self.mundane_positions = {
            # House placement, decimal longitude (out of 360*)
            "Sun" : [c_double(), c_double()],
            "Moon" : [c_double(), c_double()],
            "Mercury" : [c_double(), c_double()],
            "Venus" : [c_double(), c_double()],
            "Mars" : [c_double(), c_double()],
            "Jupiter" : [c_double(), c_double()],
            "Saturn" : [c_double(), c_double()],
            "Uranus" : [c_double(), c_double()],
            "Neptune" : [c_double(), c_double()],
            "Pluto" : [c_double(), c_double()],
            }
        self.cusps = {
            "1" : (c_double)(),
            "2" : (c_double)(),
            "3" : (c_double)(),
            "4" : (c_double)(),
            "5" : (c_double)(),
            "6" : (c_double)(),
            "7" : (c_double)(),
            "8" : (c_double)(),
            "9" : (c_double)(),
            "10" : (c_double)(),
            "11" : (c_double)(),
            "12" : (c_double)(),
            }
        self.angles = {
            "Asc" : (c_double)(),
            "MC" : (c_double)(),
            "Dsc" : (c_double)(),
            "IC" : (c_double)(),
            "EP" : (c_double)(),
            "Zen" : (c_double)(),
            "WP" : (c_double)(),
            "Ndr" : (c_double)(),
            }


# Swiss Ephemeris DLL functions, wrapped for use by Python
file_dir = os.path.dirname(os.path.abspath(__file__))
dll = windll.LoadLibrary(file_dir + '/SE/sweph/bin/swedll64.dll')

# sets the filepath of the ephemeris for the DLL functions
py_set_ephemeris_path = dll.swe_set_ephe_path
py_set_ephemeris_path.argtypes = [c_char_p]
py_set_ephemeris_path.restype = None

py_set_sidereal_mode = dll.swe_set_sid_mode
py_set_sidereal_mode.argtypes = (c_int32, c_double, c_double)
py_set_sidereal_mode.restype = None

# Convert local time, with timezone, to UTC
# Arg order: (in-year, in-month, in-day, in-hour, in-min, in-sec,
# in-timezone, out-year, out-month, out-day, out-hour, out-min, out-sec)
py_local_time_to_UTC = dll.swe_utc_time_zone
py_local_time_to_UTC.argtypes = [c_int32, c_int32, c_int32, c_int32, 
                                 c_int32, c_double, c_double, 
                                 POINTER(c_int32), POINTER(c_int32), 
                                 POINTER(c_int32), POINTER(c_int32), 
                                 POINTER(c_int32), POINTER(c_double)]
py_local_time_to_UTC.restype = None

# Calculates julian day number (in universal time) from:
# year, month, day, fractional hour, calendar flag (should be 1)
py_get_julian_day = dll.swe_julday
py_get_julian_day.argtypes = [c_int, c_int, c_int, c_double, c_int]
py_get_julian_day.restype = c_double

# Accepts the Julian Day Number in Universal Time, 
# and returns sidereal time at the Greenwich Meridian.
# Must be converted to local sidereal time for mundane calculations.
py_get_sidereal_time_UTC = dll.swe_sidtime
py_get_sidereal_time_UTC.argtypes = [c_double]
py_get_sidereal_time_UTC.restype = c_double

# Calculates planetary positions from:
# Julian Day (universal time), body #, zodiacal flag, 
# pointer to array of 6 doubles to write in, string to write errors to)
py_calculate_planets_UT = dll.swe_calc_ut
py_calculate_planets_UT.argtypes = [c_double, c_int, c_int32, 
                                    POINTER(c_double), c_char_p]
py_calculate_planets_UT.restype = None

# Calculates ayanamsa from: Julian Day in UT, sidereal flag, 
# pointer to double to write ayanamsa in, pointer to error string.
# Returns either the ephemeris flag (a positive int), or ERR (-1)
py_get_ayanamsa_UT = dll.swe_get_ayanamsa_ex_ut
py_get_ayanamsa_UT.argtypes = [c_double, c_int32, 
                               POINTER(c_double), c_char_p]
py_get_ayanamsa_UT.restype = c_int32


# Calculates Ascendant, MC, and house cusps.
# The arguments are Julian Day in UT, the ephemeris flag (sidereal),
# Geolatitude, geolongitude, int house system ('C'), double array cusps, double array AscMc (etc)
# swe_houses_ex() returns an int, but the documention doesn't specify what it is, so we retype to None
py_calculate_houses = dll.swe_houses_ex
py_calculate_houses.argtypes = [c_double, c_int32, c_double, c_double, 
                                c_int, POINTER(c_double), POINTER(c_double)]
py_calculate_houses.restype = None


# Python-level functions

def get_sign(longitude):

    zodiac = ["Ari", "Tau", "Gem", "Can", "Leo", "Vir", 
              "Lib", "Sco", "Sag", "Cap", "Aqu", "Pis"]
    key = int(longitude/30)
    return zodiac[key]

def parse_aspect(pname1, plong1, pname2, plong2):
    """Checks two planets' positions for any valid aspects, 
    returning that aspect info as a string if applicable"""

    if pname1 == pname2:
        return None
    
    orb = 0
    tier = None
    aspect_type = ""

    def get_orb(longitude_one, longitude_two, lowbound, highbound):
        """Gets the orb, or distance from exact, of an aspect"""

        aspect = math.fabs(longitude_one - longitude_two)
        aspect360 = math.fabs(aspect - 360)
        aspect_average = (lowbound + highbound) / 2

        if aspect >= lowbound and aspect <= highbound:
            if lowbound != 0:
                return math.fabs(aspect - aspect_average)
            else:
                return aspect
        elif aspect360 >= lowbound and aspect360 <= highbound:
            if lowbound != 0:
                return math.fabs(aspect360 - aspect_average)
            else:
                return aspect360

    def get_priority(pname1, pname2, _tier, atype, relative_orb):
        """Function to calculate weight for aspects; to be utilized in a later release"""
        
        priority = 0
        lights = ["Sun", "Moon"]
        outer_planets = ["Uranus", "Neptune", "Pluto"]
        hard_aspects = ["Cnj", "Opp", "Sqr", "Sms", "Ssq"]
        
        if pname1 in lights:
            priority += 1
        if pname2 in lights:
            priority += 1

        priority += (3 - _tier)

        if atype in hard_aspects:
            priority += 0.10
        if pname1 in outer_planets and pname2 in outer_planets:
            priority -= 0.5
        priority = priority - relative_orb
        return priority

    # Case: conjunction
    if get_orb(plong1, plong2, 0, 10) != None:
        aspect_type = "Cnj"
        orb = get_orb(plong1, plong2, 0, 10)
        if orb <= 1:
            tier = 0
        if orb <= 4:
            tier = 1
        elif orb <= 7:
            tier = 2
        elif orb <= 10:
            tier = 3
        priority_placeholder = get_priority(pname1, pname2, tier, aspect_type, (orb/5))
    else:
        pass
    
    # Case: opposition
    if get_orb(plong1, plong2, 170, 190) != None:
        aspect_type = "Opp"
        orb = (get_orb(plong1, plong2, 170, 190) % 180)
        if orb <= 1:
            tier = 0
        if orb <= 4:
            tier = 1
        elif orb <= 7:
            tier = 2
        elif orb <= 10:
            tier = 3
        priority_placeholder = get_priority(pname1, pname2, tier, aspect_type, (orb/5))
    else:
        pass

    # Case: square
    if get_orb(plong1, plong2, 82.5, 97.5) != None:
        aspect_type = "Sqr"
        orb = (get_orb(plong1, plong2, 82.5, 97.5) % 90)
        if orb <= 1:
            tier = 0
        if orb <= 3:
            tier = 1
        elif orb <= 6:
            tier = 2
        elif orb <= 7.5:
            tier = 3
        priority_placeholder = get_priority(pname1, pname2, tier, aspect_type, (orb/4.8))
    else:
        pass

    # Case: trine
    if get_orb(plong1, plong2, 115, 125) != None:
        aspect_type = "Tri"
        orb = (get_orb(plong1, plong2, 115, 125) % 120)
        if orb <= 1:
            tier = 0
        if orb <= 3:
            tier = 1
        elif orb <= 5:
            tier = 2
        else: 
            tier = None
        priority_placeholder = get_priority(pname1, pname2, tier, aspect_type, (orb/4.5))
    else:
        pass

    # Case: sextile
    if  get_orb(plong1, plong2, 55, 65) != None: 
        aspect_type = "Sxt"
        orb = (get_orb(plong1, plong2, 55, 65) % 60)
        if orb <= 1:
            tier = 0
        if orb <= 3:
            tier = 1
        elif orb <= 5:
            tier = 2
        else: 
            tier = None
        priority_placeholder = get_priority(pname1, pname2, tier, aspect_type, (orb/4.5))
    else:
        pass

    # Case: semisquare
    if get_orb(plong1, plong2, 43, 47) != None: 
        aspect_type = "Sms"
        orb = (get_orb(plong1, plong2, 43, 47) % 45)
        if orb <= 1:
            tier = 0
        if orb <= 2:
            tier = 1
        else: 
            tier = None
        priority_placeholder = get_priority(pname1, pname2, tier, aspect_type, (orb/4))
    else:
        pass

    # Case: sesquisquare
    if get_orb(plong1, plong2, 133, 137) != None: 
        aspect_type = "Ssq"
        orb = (get_orb(plong1, plong2, 133, 137) % 135)
        if orb <= 1:
            tier = 0
        if orb <= 2:
            tier = 1
        else: 
            tier = None
        priority_placeholder = get_priority(pname1, pname2, tier, aspect_type, (orb/4))
    else:
        pass

    # Case: quincunx
    if get_orb(plong1, plong2, 71, 73) != None: 
        aspect_type = "Qnx"
        orb = (get_orb(plong1, plong2, 71, 73) % 72)
        if orb <= 1:
            tier = 1
        else: 
            tier = None
        priority_placeholder = get_priority(pname1, pname2, tier, aspect_type, (orb/2))
    else:
        pass

    if orb != 0 and tier is not None:
        returnvalue = []
        orb_deg = math.trunc(orb)
        orb_min = math.trunc((orb - (math.trunc(orb))) * 60)
        returnvalue += pname1 + " " + aspect_type + " " + pname2 + " " + str(orb_deg) + "* " + str(orb_min) + "'"
        aspect_return = "".join(returnvalue)
        return (priority_placeholder, aspect_return)
    else:
       return None                       

def get_LST(year, month, day, decimalhour, timezone, decimal_longitude):
    """Calculates local sidereal time for date, time, location of event"""

    # Julian Day number for midnight
    julian_day_0_GMT = py_get_julian_day(year, month, day, 0, 1)

    # Local time converted to UT/GMT
    universal_time = (decimalhour - timezone)

    # Sidereal time for the Julian day number, at midnight
    sidereal_time_0_GMT = (julian_day_0_GMT - 2451545.0) / 36525.0

    # "LST" for the Greenwich Meridian, not local yet
    greenwich_sidereal_time = (6.697374558 
                               + (2400.051336 * sidereal_time_0_GMT) 
                               + (0.000024862 
                               * (math.pow(sidereal_time_0_GMT, 2))) 
                               + (universal_time * 1.0027379093))


    local_sidereal_time = ((greenwich_sidereal_time 
                            + (decimal_longitude / 15)) % 24)

    if local_sidereal_time < 0:
        local_sidereal_time += 24

    return local_sidereal_time

def calculate_mundane_positions(planet_latitude, planet_longitude, 
                                 LST, obliquity, decimal_svp, 
                                 decimal_geolat, planet_pvl):
    """ Calculate a planet's mundane angularity (house position in prime vertical longitude) """

#region calculations
    ramc = LST * 15

    calc_ax = (math.cos(math.radians(planet_longitude 
                    + (360 - (330 + decimal_svp)))))

    precessed_d = (math.degrees(math.asin
                    (math.sin(math.radians(planet_latitude)) 
                    * math.cos(math.radians(obliquity)) 
                    + math.cos(math.radians(planet_latitude)) 
                    * math.sin(math.radians(obliquity)) 
                    * math.sin(math.radians(planet_longitude 
                    + (360 - (330 + decimal_svp)))))))

    calc_ay = (math.sin(math.radians((planet_longitude 
                    + (360 - (330 + decimal_svp))))) 
                    * math.cos(math.radians(obliquity)) 
                    - math.tan(math.radians(planet_latitude)) 
                    * math.sin(math.radians(obliquity)))

    calc_ayx_deg = math.degrees(math.atan(calc_ay / calc_ax))

    precessed_a = None
    if (calc_ax < 0):
        precessed_a = calc_ayx_deg + 180
    else: 
        if(calc_ay < 0):
            precessed_a = calc_ayx_deg + 360
        else: precessed_a = calc_ayx_deg

    ha_deg = ramc - precessed_a

    calc_cz = (math.degrees(math.atan(1 
             / (math.cos (math.radians(decimal_geolat)) 
             / math.tan(math.radians(ha_deg)) 
             + math.sin(math.radians(decimal_geolat)) 
             * math.tan(math.radians(precessed_d)) 
             / math.sin(math.radians(ha_deg))))))

    calc_cx = (math.cos(math.radians(decimal_geolat)) 
               * math.cos(math.radians(ha_deg)) 
               + math.sin(math.radians(decimal_geolat)) 
               * math.tan(math.radians(precessed_d)))

    campanus_longitude = ()
    if (calc_cx < 0):
        campanus_longitude = 90 - calc_cz
    else:
       campanus_longitude = 270 - calc_cz

#endregion

    planet_pvl[0] = (int(campanus_longitude/30) + 1)
    planet_pvl[1] = campanus_longitude
    return None

def calculate_houses(julian_day_number, geolatitude, geolongitude,
                     classname):
    """Calculates house cusps and ecliptical 
    longitudes of angles in the Campanus system"""

    # Swiss Ephemeris documentation specifies 8 doubles for houses, 13 for cusps
    cusp_array = (c_double * 13)()
    house_array = (c_double * 8)()
    py_calculate_houses(julian_day_number, SIDEREALMODE,
                        geolatitude, geolongitude, CAMPANUS,
                        cusp_array, house_array)
    classname.angles["Asc"] = house_array[0]
    classname.angles["MC"] = house_array[1]
    classname.angles["Dsc"] = (house_array[0] + 180) % 360
    classname.angles["IC"] = (house_array[1] + 180) % 360
    classname.angles["EP"] = (classname.angles["MC"] + 90) % 360
    classname.angles["Zen"] = (classname.angles["Dsc"] + 90) % 360
    classname.angles["WP"] = (classname.angles["IC"] + 90) % 360
    classname.angles["Ndr"] = (classname.angles["Asc"] + 90) % 360    

    return None

def calculate_foreground_planets(classname):
    """Calculates foreground and background planets, and returns them as a tuple"""

    foregroundlist = []
    backgroundlist = []
    
    # Possible house locations that would make a planet foreground 
    primary_angles = {
        "Asc": (12, 1), 
        "MC" : (9, 10), 
        "Dsc" : (6, 7), 
        "IC" : (3, 4)
        }
    secondary_angles = ["EP", "Zen", "WP", "Ndr"]
    background_angles = [(2, 3), (5, 6), (8, 9), (11, 12)]

    # Planet proximity to primary angles, measured in prime vertical longitude
    for angle in primary_angles.keys():
        for planet in classname.mundane_positions.keys():
            house = classname.mundane_positions[planet][0]
            longitude = classname.mundane_positions[planet][1] % 30
            if ((house == primary_angles[angle][0] and longitude >= 20) 
                or (house == primary_angles[angle][1] and longitude <= 10)):

                if longitude >= 20:
                    orb = math.fabs(longitude - 30)
                else:
                    orb = longitude
                returnvalue = []
                orb_deg = math.trunc(orb)
                orb_min = math.trunc((orb - (math.trunc(orb))) * 60)
                returnvalue += planet + " Cnj " + angle + " " + str(orb_deg) + "* " + str(orb_min) + "'"
                foreground_planet = "".join(returnvalue)
                foregroundlist.append(foreground_planet)    

    # Planet proximity to background cusps, measured in prime vertical longitude
    for angle in background_angles:
        for planet in classname.mundane_positions.keys():
            house = classname.mundane_positions[planet][0]
            longitude = classname.mundane_positions[planet][1] % 30
            if (house == angle[0] and longitude >= 20) or (house == angle[1] and longitude <= 10):
                if longitude >= 20:
                    orb = math.fabs(longitude - 30)
                else:
                    orb = longitude
                returnvalue = []
                orb_deg = math.trunc(orb)
                orb_min = math.trunc((orb - (math.trunc(orb))) * 60)
                returnvalue += planet + " background " + str(orb_deg) + "* " + str(orb_min) + "'"
                background_planet = "".join(returnvalue)
                backgroundlist.append(background_planet)

    # Planet proximity to secondary angles, measured in ecliptical longitude
    for angle in secondary_angles:
        for key in PLANETLIST:
            ecliptical_longitude = classname.planet_dictionary[key][0]
            point = classname.angles[angle]

            # So we don't have to deal with longitudes near 360* as such
            if ecliptical_longitude >= 355:
                ecliptical_longitude -= 360
            if point >= 355:
                point -= 360
        
            if ecliptical_longitude >= point - 3 and ecliptical_longitude <= point + 3:
                
                returnvalue = []
                orb = math.fabs(ecliptical_longitude - point)
                orb_deg = math.trunc(orb)
                orb_min = math.trunc((orb - (math.trunc(orb))) * 60)
                returnvalue += key + " Cnj " + angle + " " + str(orb_deg) + "* " + str(orb_min) + "'"
                foreground_planet = "".join(returnvalue)
                foregroundlist.append(foreground_planet)
 
    return (foregroundlist, backgroundlist)

def calculate_natal_data(year, month, day, hour, 
                         minute, second, utc_offset, decimal_longitude, 
                         decimal_latitude, classname):
    """ Populates class instance with birth info, planetary coordinates """

    classname.datetimelocation["Year"] = year
    classname.datetimelocation["Month"] = month
    classname.datetimelocation["Day"] = day
    classname.datetimelocation["Hour"] = hour
    classname.datetimelocation["Min"] = minute
    classname.datetimelocation["Sec"] = second
    classname.datetimelocation["Longitude"] = decimal_longitude
    classname.datetimelocation["Latitude"] = decimal_latitude

    decimalhour_local = (((((hour * 60)      
                         + minute) * 60)      
                         + second) / 3600)
    
    LST = get_LST(year, month, day, decimalhour_local, 
                  utc_offset, decimal_longitude)

    # Conversion from local time to UTC
    (outyear, outmonth, outday, 
     outhour, outmin, outsec) = (c_int32(), c_int32(), c_int32(), 
                                 c_int32(), c_int32(), c_double())
    py_local_time_to_UTC(year, month, day, hour, minute, second, 
                         utc_offset, outyear, outmonth, outday, outhour, 
                         outmin, outsec)

    decimalhour_UTC = (((((outhour.value * 60)      
                         + outmin.value) * 60)      
                         + outsec.value) / 3600)    

    # Using the new UTC time to get a correct Julian Day Number
    time_julian_day = py_get_julian_day(outyear, outmonth, outday, 
                                        decimalhour_UTC, 1)

    errorstring = create_string_buffer(126)

    SVP = c_double()
    ayanamsa_return = py_get_ayanamsa_UT(time_julian_day, 
                                         SIDEREALMODE, SVP, errorstring)
    if (ayanamsa_return < 0):
        print("Error retrieving ayanamsa")

    # Difference in tropical and sidereal positions desired,
    # not the longitude of one point in the other zodiac 
    SVP = (30 - SVP.value)

    planet_number = 0
    returnarray = [(c_double * 6)() for x in range(10)]
    for key in classname.planet_dictionary.keys():

        # dll.swe_calc_ut uses ints as identifiers; 0-9 is Sun-Pluto
        if planet_number <= 9:
            py_calculate_planets_UT(time_julian_day, 
                                    planet_number, SIDEREALMODE, 
                                    returnarray[planet_number], 
                                    errorstring) 
            classname.planet_dictionary[key] = returnarray[planet_number]
            planet_number += 1

        else:

            # -1 is the special "planetary body" for calculating obliquity
            py_calculate_planets_UT(time_julian_day, -1, SIDEREALMODE, 
                                    classname.planet_dictionary["Obliquity"], 
                                    errorstring)
            break

    for key in classname.mundane_positions.keys():
        calculate_mundane_positions(classname.planet_dictionary[key][1], 
                                     classname.planet_dictionary[key][0], 
                                     LST, 
                                     classname.planet_dictionary["Obliquity"][0], 
                                     SVP, decimal_latitude, 
                                     classname.mundane_positions[key])
    
    classname.planet_dictionary["SVP"] = SVP
    classname.planet_dictionary["LST"] = LST
    classname.planet_dictionary["Longitude"] = decimal_longitude
    classname.planet_dictionary["Latitude"] = decimal_latitude

    calculate_houses(time_julian_day, decimal_latitude, decimal_longitude, classname)

    return None

def calculate_and_sort_aspects(classname):
    """Prioritizes the order of aspects based on their impact 
    in the chart, as deterined by get_priority()"""

    aspect_list = []
    aspect_priority = 0
    loopcount = 0 
    for key in PLANETLIST:
        for planetcounter in range (loopcount, 10):
            planetname = str(PLANETLIST[planetcounter])
            potential_aspect = (parse_aspect(str(key), classname.planet_dictionary[key][0], planetname,
                        classname.planet_dictionary[planetname][0]))
            if potential_aspect is not None:
                aspect_list.append(potential_aspect)
            planetcounter += 1
        loopcount += 1
    aspect_list.sort(reverse=True) 
    return aspect_list

def print_natal_data(classname):
    """Prints all relevant natal information to a .txt file"""

    natalfile = open("{}.txt".format(classname.name), "w+")
    natalfile.write("~* AstroNova v. {} *~\n".format(VERSION_NUMBER))
    natalfile.write("Natal instance: {}\n".format(classname.name))
    natalfile.write("{} {} {} {}:{}\n".format(classname.datetimelocation["Year"], 
                                            classname.datetimelocation["Month"], 
                                            classname.datetimelocation["Day"], 
                                            classname.datetimelocation["Hour"], 
                                            classname.datetimelocation["Min"]))
    natalfile.write("Long: {}   Lat: {}\n\n\n".format(classname.datetimelocation["Longitude"], 
                                                classname.datetimelocation["Latitude"]))
    natalfile.write("Sign Placements: \n\n")
    for key in PLANETLIST:
        natalfile.write("{}: {}*{}' {}".format(key, 
            int((classname.planet_dictionary[key][0] % 30)), 
            (int(round(((classname.planet_dictionary[key][0] % 30) 
            - math.floor((classname.planet_dictionary[key][0] % 30))) * 60))), 
            get_sign(classname.planet_dictionary[key][0])))
        natalfile.write("\n")
    natalfile.write("\n")

    natalfile.write("\nPrimary Angles (Eclipto): \n\n")
    for key in classname.angles.keys():
        if key == "Asc" or key == "Dsc" or key == "MC" or key == "IC":
            natalfile.write("{} {}* {}' {}".format(key, 
                                            ((math.floor(classname.angles[key])) % 30), 
                                            round((classname.angles[key] 
                                            - math.floor(classname.angles[key])) * 60), 
                                            get_sign(classname.angles[key])))
            natalfile.write("\n")
    natalfile.write("\n")

    # for debugging, uncomment the code below
    #natalfile.write("\nMundane House Placements: \n\n")
    #for key in PLANETLIST:
    #    natalfile.write("{} house position: {}, {}* {}'".format(key, classname.mundane_positions[key][0], 
    #                                                    int(classname.mundane_positions[key][1] 
    #                                                        - (int(classname.mundane_positions[key][1]/30)) * 30), 
    #                                                    (int((classname.mundane_positions[key][1] 
    #                                                    - int(classname.mundane_positions[key][1])) * 60))))
    #    natalfile.write("\n")
    
    angularity_list = calculate_foreground_planets(classname)

    natalfile.write("\nForeground Planets: \n\n")
    for placement in angularity_list[0]:
        natalfile.write(placement)
        natalfile.write("\n")

    natalfile.write("\n")

    natalfile.write("\nBackground Planets: \n\n")
    for placement in angularity_list[1]:
        natalfile.write(placement)
        natalfile.write("\n")
    natalfile.write("\n")

    natalfile.write("\nList of Aspects: \n\n")
    sorted_aspect_list = calculate_and_sort_aspects(classname)
    for priority, aspect in sorted_aspect_list:
        natalfile.write("{} -- priority: {}\n".format(aspect, priority))
    natalfile.close()
