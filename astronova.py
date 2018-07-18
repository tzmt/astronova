from astrosetup import *
geolocator = Nominatim()
timezonefinder = TimezoneFinder()

def main():

    def interface_submission(*args):
        """Primary function called submitting birth data via interface"""

        native_name = None
        location = None
        birthdate = None
        birthtime = None
        splitbirthdate = None
        ssr_year = None
        slr_start = None
        relocation = None
        radiovalue = None
        native_instance = None
        local_instance = None
        ssr_instance = None
        slr_instance = None

        try:
            native_name = str(native_name_value.get())
            native_name = native_name.strip(" ")
            if native_name == "":
                messagebox.showinfo(message="Invalid name")
                return None
        except:
            raise RuntimeError("Unable to retrieve name")

        try:
            location = str(location_value.get())
            location = location.strip()
            location = location.lower()
            if native_name == "test":
                location = "ridgewood, nj, usa"
            elif location == "":
                messagebox.showinfo(message="Invalid location")
                return None

        except:
            raise RuntimeError("Error retrieving location")

        try:
            birthdate = str(birthdate_value.get())
            birthdate = birthdate.strip()

            if native_name == "test":
                birthdate = "12/20/1989"
            elif birthdate == "":
                messagebox.showinfo(message="Invalid birthdate")
                return None

        except:
            raise RuntimeError("Error retrieving birthdate")
        
        try:
            birthtime = str(birthtime_value.get())
            birthtime = birthtime.strip()
            if native_name == "test":
                birthtime = "10:20pm"
            elif birthtime == "":
                messagebox.showinfo(message="Invalid birth time")
                return None

        except:
            raise RuntimeError("Error retrieving birth time")

        try:
            splitbirthdate = birthdate.split("/")
        except:
            raise RuntimeError("Error splitting birthdate")

        if len(splitbirthdate) == 3:
            month, day, year = (x for x in splitbirthdate)
            # the below line can probably be removed, and be added to the above assignment generator: int(x) for x in splitbirthddate
            month, day, year = int(month), int(day), int(year)
        else:
            messagebox.showinfo(message="Invalid format; use month/day/year")
            return None
        if year < 0 or year > 2100:
            messagebox.showinfo(message="Year out of range; choose year >= 0 and <= 2100")
            return None
        elif year < 30:
            year += 2000
        elif year < 100 and year >= 30:
            year += 1900
        if month < 1 or month > 12:
            messagebox.showinfo(message="Invalid month")
            return None
        if day < 1 or day > 31:
            messagebox.showinfo(message="Invalid day")
            return None

        try:
            ampm = (birthtime.strip("0123456789: ")).lower()
            birthtime = birthtime.strip("AMPamp ")
            if ":" not in birthtime:
                if len(birthtime) == 2:
                    hour = birthtime
                    min = 0
                elif len(birthtime) == 3:
                    hour = birthtime[0]
                    min = birthtime[1:3]
                elif len(birthtime) == 4:
                    hour = int(birthtime[0:2])
                    min = birthtime[2:4]
            else: 
                hour, min = birthtime.split(":")

            hour, min = int(hour), int(min)

        except:
            raise RuntimeError("Error splitting birthtime")
        if hour < 0 or hour > 23:
            messagebox.showinfo(message="Invalid hour")
            return None
        if min < 0 or min > 59:
            messagebox.showinfo(message="Invalid minute")
            return None

        if "a" in str(ampm) and "p" in str(ampm):
            messagebox.showinfo(message="Invalid am/pm information")
        elif "p" in str(ampm) and hour < 12:
            hour += 12
        elif "a" in str(ampm) and hour == 12:
            hour = 0
        elif "a" in str(ampm) and hour > 12:
            messagebox.showinfo(message="Invalid birth time")
            return None

        try:
            relocation = str(relocation_value.get())
            relocation = relocation.strip()
            if native_name == "test":
                relocation = "hackensack, nj, usa"
            elif relocation == "":
                relocation = None
        except:
            raise RuntimeError("Error retrieving local residence")

        try:
            radiovalue = (radiobutton_value.get()).strip()
            if radiovalue == "SSR":
                try:
                    ssr_year = (solunar_value.get()).strip()
                    if ssr_year != "":
                        if not ssr_year.isnumeric():                     
                            messagebox.showinfo(message="Invalid SSR year")
                    else:
                        ssr_year = None
                except:
                    raise RuntimeError("Error retrieving SSR year")

            elif radiovalue == "SLR":
                try:
                    slr_start = (str(solunar_value.get())).strip()
                    if slr_start == "":
                        slr_start = None
                    else:
                        try:
                            slr_start_split = slr_start.split("/")
                        except:
                            raise RuntimeError("Error splitting SLR start date")
                except:
                    raise RuntimeError("Error retrieving SLR start date")


            elif radiovalue == "DSLR":
                try:
                    slr_start = (str(solunar_value.get())).strip()
                    if slr_start == "":
                        slr_start = None
                    else:
                        try:
                            slr_start_split = slr_start.split("/")
                        except:
                            raise RuntimeError("Error splitting DSLR date")
                except:
                    raise RuntimeError("Error retrieving DSLR date")

            elif radiovalue == "Transits":
                try:
                    slr_start = (str(solunar_value.get())).strip()
                    if slr_start == "":
                        slr_start = None
                    else:
                        try:
                            slr_start_split = slr_start.split("/")
                        except:
                            raise RuntimeError("Error splitting transit date")
                except:
                    raise RuntimeError("Error retrieving transit date")
            
            else:
                ssr_year, slr_start = None, None
        except:
            raise RuntimeError("Error determining chart type!")

        if slr_start is not None:
            if len(slr_start_split) == 3:
                startmonth, startday, startyear = (int(x) for x in slr_start_split)

            else:
                messagebox.showinfo(message="Invalid format; use month/day/year")
                return None
            if startyear < 0 or startyear > 2100:
                messagebox.showinfo(message="Year out of range; choose year >= 0 and <= 2100")
                return None
            elif startyear < 30:
                startyear += 2000
            elif startyear < 100 and startyear >= 30:
                startyear += 1900
            if startmonth < 1 or startmonth > 12:
                messagebox.showinfo(message="Invalid month")
                return None
            if startday < 1 or startday > 31:
                messagebox.showinfo(message="Invalid day")
                return None

        native_instance = Natal(native_name)
        try:
            natal_location = geolocator.geocode(location)
        except:
            messagebox.showinfo(message="Natal location lookup failure: Have you tried asking nicely?")
            return None
        if natal_location is None:
            messagebox.showinfo(message="Unable to locate region; try different name, format, or city")
            return None
        natal_latitude = natal_location.latitude
        natal_longitude = natal_location.longitude
        if natal_latitude is None or natal_longitude is None:
            messagebox.showinfo(message="Natal latitude or longitude are unavailable; try a different region")
            return None
        else:
            natal_timezone = timezonefinder.timezone_at(lng=natal_longitude, lat=natal_latitude)

            native_instance.datetime = pendulum.create(int(year), int(month), int(day), int(hour), int(min), 0, 0, natal_timezone)
            natal_utc_offset = native_instance.datetime.offset_hours

        local_instance = None
        if relocation is not None:
            try:
                local_residence = geolocator.geocode(relocation)
            except:
                messagebox.showinfo(message="Lookup failure for current residence; have you tried turning it off and on again?")
                return None
            if local_residence is None:
                messagebox.showinfo(message="Unable to locate lat/long for current residence; try different name, format, or city")
                return None
            local_latitude = local_residence.latitude
            local_longitude = local_residence.longitude
            if local_latitude is None or local_longitude is None:
                messagebox.showinfo(message="Local latitude or longitude are unavailable; try a different region")
                return None
            local_timezone = timezonefinder.timezone_at(lng=local_longitude, lat=local_latitude)
            print("Relocation successful! Using {} timezone".format(local_timezone))
            print(local_residence)
            local_tz_pendulum = pendulum.timezone(local_timezone)
            relocated_datetime = local_tz_pendulum.convert(native_instance.datetime)
            native_instance.datetime = relocated_datetime
            print("Adjusted birth time: {}".format(native_instance.datetime.to_day_datetime_string()))
            local_instance = native_instance
            local_instance.location["Longitude"] = local_longitude 
            local_instance.location["Latitude"] = local_latitude
            print(local_residence.longitude)
            print(local_residence.latitude)
            print(local_longitude)
            print(local_latitude)

        native_instance.location["Longitude"] = local_longitude if relocation != None else natal_longitude 
        native_instance.location["Latitude"] = local_latitude if relocation != None else natal_latitude

        # Sub-minute accuracy is not required for calculation of natal data
        calculate_natal_data(native_instance if local_instance == None else local_instance)
        print_chart_data(native_instance, "Natal")
        messagebox.showinfo(message="Calculation successful! See {}.txt in the AstroNova program folder".format(native_name))

        if radiovalue == "SSR" and ssr_year is not None:
            ssr_year = int(ssr_year)
            ssr_instance = Natal("{} SSR".format(native_instance.name))
            ssr_instance.location["Longitude"] = native_instance.location["Longitude"] if relocation == None else local_longitude
            ssr_instance.location["Latitude"] = native_instance.location["Latitude"] if relocation == None else local_latitude
            ssr_instance.datetime = native_instance.datetime.add(years=(ssr_year - native_instance.datetime.year))
            calculate_ssr_chart(native_instance if local_instance == None else local_instance, ssr_instance)
            print_full_solunar_return(native_instance if local_instance == None else local_instance, ssr_instance, radiovalue)


        if (radiovalue == "SLR" or radiovalue == "DSLR") and slr_start is not None:
            startdate = pendulum.create(startyear, startmonth, startday, 0, 0, 0, 0, natal_timezone)
            enddate = startdate.add(days=14)
            slr_instance = Natal("{} {}".format(native_instance.name, radiovalue))
            slr_instance.location["Longitude"] = native_instance.location["Longitude"] if relocation == None else local_longitude
            slr_instance.location["Latitude"] = native_instance.location["Latitude"] if relocation == None else local_latitude
            slr_instance.datetime = native_instance.datetime
            calculate_slr_chart(native_instance if local_instance == None else local_instance, slr_instance, startdate, enddate, radiovalue)
            print_full_solunar_return(native_instance if local_instance == None else local_instance, slr_instance, radiovalue)

        if radiovalue == "Transits":

            ssr_year = startyear if native_instance.datetime.month < startmonth else startyear - 1
            ssr_instance = Natal("{} SSR".format(native_instance.name))
            ssr_instance.location["Longitude"] = native_instance.location["Longitude"]
            ssr_instance.location["Latitude"] = native_instance.location["Latitude"]
            ssr_instance.datetime = native_instance.datetime.add(years=(ssr_year - native_instance.datetime.year))
            calculate_ssr_chart(native_instance if local_instance == None else local_instance, ssr_instance)

            transit_date = pendulum.create(startyear, startmonth, startday, 0, 0, 0, 0, native_instance.datetime.timezone)
            transit_instance = Natal("{} Transits".format(native_instance.name))
            transit_instance.location["Longitude"] = local_instance.location["Longitude"] if relocation != None else native_instance.location["Longitude"]
            transit_instance.location["Latitude"] = local_instance.location["Latitude"] if relocation != None else native_instance.location["Latitude"]
            transit_instance.datetime = transit_date

            # None of this makes any sense at all
            #prog_jul_day = get_sp_julian_day(native_instance, transit_date)
            #prog_instance = Natal("Prog")
            #calculate_natal_data(prog_instance)

            print_active_transits(native_instance, local_instance, ssr_instance, transit_instance, transit_date)


        return None


    epath = os.path.dirname(os.path.abspath(__file__))                              
    epath = os.path.join(epath + "\\SE\\sweph\\ephemeris\\")
    epath = epath.encode('utf-8')
    epointer = c_char_p(epath)
    py_set_ephemeris_path(epointer)

    py_set_sidereal_mode (0, 0, 0)
    
    # Create the interface
    root = Tk()
    root.title("AstroNova v{}".format(VERSION_NUMBER))
    mainframe = ttk.Frame(root, borderwidth=5, padding=(40, 40, 40, 40))
    mainframe.grid(column=0, row=0, sticky=(N, W, E, S))
    mainframe.columnconfigure(0, weight=1)
    mainframe.rowconfigure(0, weight=1)

    native_name_value = StringVar()
    name_label = ttk.Label(mainframe, text="Name ")
    name_label.grid(column=1, columnspan=2, sticky=E, row=0, padx=2, pady=2)

    name_entry = ttk.Entry(mainframe, width=24, textvariable=native_name_value)
    name_entry.grid(column=3, row=0, padx=2, pady=2)

    location_value = StringVar()
    location_label = ttk.Label(mainframe, text="City, State/Country ")
    location_label.grid(column=1, columnspan=2, sticky=E, row=2, padx=2, pady=2)

    location_entry = ttk.Entry(mainframe, width=24, textvariable=location_value)
    location_entry.grid(column=3, row=2, padx=2, pady=2)

    birthdate_value = StringVar()
    birthdate_label = ttk.Label(mainframe, text="Date (mm/dd/yyyy)")
    birthdate_label.grid(column=1, columnspan=2, sticky=E, row=4, padx=2, pady=2)

    birthdate_entry = ttk.Entry(mainframe, width=24, textvariable=birthdate_value)
    birthdate_entry.grid(column=3, row=4, padx=2, pady=2)

    birthtime_value = StringVar()
    birthtime_label = ttk.Label(mainframe, text="Time (hh:mm am/pm)")
    birthtime_label.grid(column=1, columnspan=2, sticky=E, row=6, padx=2, pady=2)

    birthtime_entry = ttk.Entry(mainframe, width=24, textvariable=birthtime_value)
    birthtime_entry.grid(column=3, row=6, padx=2, pady=2)

    solunar_value = StringVar()
    solunar_label = ttk.Label(mainframe, text="SSR Year / SLR Date")
    solunar_label.grid(column=2, row=9, padx=2, pady=2)

    solunar_entry = ttk.Entry(mainframe, width=24, textvariable=solunar_value)
    solunar_entry.grid(column=3, row=9, padx=2, pady=2)

    relocation_value = StringVar()
    relocation_label = ttk.Label(mainframe, text="Local Residence")
    relocation_label.grid(column=2, row=10, padx=2, pady=2)

    relocation_entry = ttk.Entry(mainframe, width=24, textvariable=relocation_value)
    relocation_entry.grid(column=3, row=10, padx=2, pady=2)

    radiobutton_value = StringVar(None, "Natal")

    Radiobutton(mainframe, text="Natal", variable=radiobutton_value, value="Natal").grid(column=1, row=7)
    Radiobutton(mainframe, text="SSR", variable=radiobutton_value, value="SSR").grid(column=2, row=7)
    Radiobutton(mainframe, text="SLR", variable=radiobutton_value, value="SLR").grid(column=3, row=7)
    Radiobutton(mainframe, text="DSLR", variable=radiobutton_value, value="DSLR").grid(column=4, row=7)
    Radiobutton(mainframe, text="Transits", variable=radiobutton_value, value="Transits").grid(column=5, row=7)

    submit_button = ttk.Button(mainframe, text="Calculate", command=interface_submission)
    submit_button.grid(column=2, row=12, columnspan=2, padx=10, pady=10)

    name_entry.focus()
    root.mainloop()

    # Free memory allocated by the DLL for raw calculations
    dll.swe_close()

main()