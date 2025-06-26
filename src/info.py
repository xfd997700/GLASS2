import pyfiglet
KG_NAME = "GLASS"
KG_VERSION = "2.0"
LOGO_FONT = "bulbhead"
VERSION_FONT = "digital"

# prints
ART_LOGO = pyfiglet.figlet_format(KG_NAME, font=LOGO_FONT)
# ART_VERSION = pyfiglet.figlet_format("v." + KG_VERSION, font=VERSION_FONT)
ART_VERSION = "ver." + KG_VERSION
PRINT_INFO = """
==============================================================
                                                      ( (
                     WELCOME!                          ) )
      Start data downloading and processing.        ........
         This process may take a while...           |      |]
     Please make yourself a coffee and chill!       \      /
                                                     `----'
==============================================================
"""
