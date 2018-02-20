from astrosetup import *
geolocator = Nominatim()
timezonefinder = TimezoneFinder()

def main():

    def interface_submission(*args):
        """Called upon submitting birth data via interface"""

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
            if location == "":
                messagebox.showinfo(message="Invalid location")
                return None
        except:
            raise RuntimeError("Error retrieving location")

        try:
            birthdate = str(birthdate_value.get())
            birthdate = birthdate.strip()
            if birthdate == "":
                messagebox.showinfo(message="Invalid birthdate")
                return None
        except:
            raise RuntimeError("Error retrieving birthdate")
        
        try:
            birthtime = str(birthtime_value.get())
            birthtime = birthtime.strip()
            if birthtime == "":
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
            month, day, year = int(month), int(day), int(year)
        else:
            messagebox.showinfo(message="Invalid format; use month/day/year")
            return None
        if year < 0 or year > 2100:
            messagebox.showinfo(message="Year out of range; choose year >= 0 and <= 2100")
            return None
        if month < 1 or month > 12:
            messagebox.showinfo(message="Invalid month")
            return None
        if day < 1 or day > 31:
            messagebox.showinfo(message="Invalid day")
            return None

        try:
            ampm = (birthtime.strip("0123456789: ")).lower()
            birthtime = birthtime.strip("AMPamp ")
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
        if str(ampm) == "pm" and hour < 12:
            hour = int(hour) + 12
        elif str(ampm) == "am" and hour == 12:
            hour = 0
        elif str(ampm) == "am" and hour > 12:
            messagebox.showinfo(message="Invalid birth time")
            return None

        native_instance = Natal(native_name)
        try:
            natal_location = geolocator.geocode(location)
        except:
            messagebox.showinfo(message="lookup failure: Unable to locate region; try different name, format, or city")
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
            natal_utc_offset = pendulum.from_timestamp(0, natal_timezone).offset_hours

        if (pendulum.create(int(year), int(month), int(day), tz=natal_timezone).is_dst) == True:
            natal_utc_offset += 1

        # Sub-minute accuracy is not required for calculation of natal data
        sec = 0
        calculate_natal_data(int(year), int(month), int(day), int(hour), 
                                int(min), sec, natal_utc_offset, natal_location.longitude, 
                                natal_location.latitude, native_instance)
        print_natal_data(native_instance)
        messagebox.showinfo(message="Calculation successful! See {}.txt in the AstroNova program folder".format(native_name))
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
    name_label.grid(column=2, row=0, padx=2, pady=2)

    name_entry = ttk.Entry(mainframe, width=24, textvariable=native_name_value)
    name_entry.grid(column=3, row=0, padx=2, pady=2)

    location_value = StringVar()
    location_label = ttk.Label(mainframe, text="City, State/Country ")
    location_label.grid(column=2, row=2, padx=2, pady=2)

    location_entry = ttk.Entry(mainframe, width=24, textvariable=location_value)
    location_entry.grid(column=3, row=2, padx=2, pady=2)

    birthdate_value = StringVar()
    birthdate_label = ttk.Label(mainframe, text="Date (mm/dd/yyyy)")
    birthdate_label.grid(column=2, row=4, padx=2, pady=2)

    birthdate_entry = ttk.Entry(mainframe, width=24, textvariable=birthdate_value)
    birthdate_entry.grid(column=3, row=4, padx=2, pady=2)

    birthtime_value = StringVar()
    birthtime_label = ttk.Label(mainframe, text="Time (hh:mm am/pm)")
    birthtime_label.grid(column=2, row=6, padx=2, pady=2)

    birthtime_entry = ttk.Entry(mainframe, width=24, textvariable=birthtime_value)
    birthtime_entry.grid(column=3, row=6, padx=2, pady=2)

    submit_button = ttk.Button(mainframe, text="Calculate Natal", command=interface_submission)
    submit_button.grid(column=2, row=12, columnspan=2, padx=10, pady=10)

    name_entry.focus()
    root.mainloop()
    
    # Free memory allocated by the DLL for raw calculations
    dll.swe_close()

main()