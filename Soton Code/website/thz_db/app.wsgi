import sys, os, bottle

sys.path = ['/var/www/thz_db'] + sys.path
os.chdir(os.path.dirname(__file__))

import main

application = bottle.default_app()
