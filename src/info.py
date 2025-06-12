import pyfiglet
KG_NAME = "GLASS"
KG_VERSION = "2.0"
LOGO_FONT = "bulbhead"
VERSION_FONT = "digital"

# prints
ART_LOGO = pyfiglet.figlet_format(KG_NAME, font=LOGO_FONT)
# ART_VERSION = pyfiglet.figlet_format("v." + KG_VERSION, font=VERSION_FONT)
ART_VERSION = "v." + KG_VERSION
PRINT_DOWNLOAD = """
==============================================================
                                                      ( (
                     WELCOME!                          ) )
    Start data downloading and pre-processing.      ........
           This process may take hours...           |      |]
     Please make yourself a coffee and chill!       \      /
                                                     `----'
==============================================================
"""

PRINT_COMPILE = """
==============================================================
                                                      ( (
                    WELCOME!                           ) )
            Start compiling metadata.               ........
          This process may take minutes...          |      |]
          Please wait a while and chill!            \      /
                                                     `----'
==============================================================
"""