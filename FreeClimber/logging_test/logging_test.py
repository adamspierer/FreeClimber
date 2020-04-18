#!/usr/bin/env python

import logging

logging.basicConfig(filename='/Users/aspierer/Manuscripts/climbing/FreeClimber/test.log',level=logging.DEBUG)
logging.debug('This message should go to the log file')
logging.info('So should this')
logging.warning('And this, too')