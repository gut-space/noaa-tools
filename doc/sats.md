NOAA 18 - AVHRR-3 sensor
- vertical half angle: 54.3 deg
- horizontal half angle: 2 deg
- 

Name 		: AVHRR/3
Role 		: OpticalImager
Health 		: Unknown
Operations 	: Advanced Very High Resolution Radiometer/3 (AVHRR/3). Refer to AVHRR/3 on NOAA 15.
Manufacturer 	: ITT Aerospace
Optical Operation # 1 
SPECTRUM CHARACTERISTICS:
# of spectrum bands: 6
Upper Spectrum Type      	: Infrared
Lower Spectrum Type      	: Visible
Spectrum Band      	:     Lower Wavelength 	 : 580 [Nanometers]    Upper Wavelength 	 : 680 [Nanometers]
Spectrum Band      	:     Lower Wavelength 	 : 725 [Nanometers]    Upper Wavelength 	 : 1000 [Nanometers]
Spectrum Band      	:     Lower Wavelength 	 : 1580 [Nanometers]   Upper Wavelength 	 : 1640 [Nanometers]
Spectrum Band      	:     Lower Wavelength 	 : 3550 [Nanometers]   Upper Wavelength 	 : 3930 [Nanometers]
Spectrum Band      	:     Lower Wavelength 	 : 1030 [Nanometers]   Upper Wavelength 	 : 1130 [Nanometers]
Spectrum Band      	:     Lower Wavelength 	 : 1150 [Nanometers]   Upper Wavelength 	 : 1250 [Nanometers]
  RESOLUTIONS:
      Resolution 1100 [Meters]
Pattern type	: Rectangular


Processing observation 2135:

noaa-apt:
cargo build && target/debug/noaa-apt --tle tle.txt --start-time "2020-08-09T05:54:02+00:00" --map yes 2135.wav

noaa-tools:
python -m noaatools.noaatools --file data/2135.png --tle tle.txt --aos "2020-08-09 05:54:02" --los "2020-08-09 06:05:40"

