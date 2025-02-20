from bottle import route, run, template, static_file, post, request, get, response
from json import dumps
from datetime import date
import extraction
import peewee

# This is the database part which is not implemented, first we define the database
# to be opend as pwdb, the data base name here is a local sql database named 'test'
# which is writtable by user 'THz' who has the passwd 'Photo-Dember'.
pwdb = peewee.MySQLDatabase('test', user = 'THz', passwd = 'Photo-Dember')

# This class defines the structure of the table with the same name as it
class THzData(peewee.Model):
    Title = peewee.TextField()
    Author = peewee.TextField()
    Description = peewee.TextField()
    Thickness = peewee.FloatField()
    xRef = peewee.TextField()
    yRef = peewee.TextField()
    xSam = peewee.TextField()
    ySam = peewee.TextField()
    nReal = peewee.TextField()
    nImag = peewee.TextField()
    Freqency = peewee.TextField()
    Date = peewee.DateField()

    class Meta:
        database = pwdb

#These are just so we can get files from subdirectories.
@route('/flot/:filename')
def flot(filename):
    return static_file(filename, root='./flot/')

@route('/jquery/:filename')
def jquery(filename):
    return static_file(filename, root='./jquery/')


@route('/js/:filename')
def js(filename):
    return static_file(filename, root='./js/')

@route('/css/:filename')
def css(filename):
    return static_file(filename, root='./css/')

@route('/css/images/:filename')
def images(filename):
    return static_file(filename, root='./css/images/')

# path route returns tabs.tpl which sets up the tab system
@route('/')
def index():
    return template('tabs')

# Post to add data to database. This is a form post with the fields as described.
@post('/addToDB')
def do_addToDB():
    newData = THzData()
    newData.Title = request.forms.get('title')
    newData.Author = request.forms.get('author')
    newData.Description = request.forms.get('description')
    newData.xRef = request.forms.get('xRVal')
    newData.yRef = request.forms.get('yRVal')
    newData.xSam = request.forms.get('xSVal')
    newData.ySam = request.forms.get('ySVal')
    thickness = request.forms.get('thick')
    newData.nReal = request.forms.get('nReal')
    newData.nImag = request.forms.get('nImag')
    newData.Freqency = request.forms.get('freq')
    newData.Date = date.today()    
    newData.Thickness = float(thickness)
    newData.save()
    return   

# Upload data tab, this tab is for uploading and displaying the time domain data
@get('/upData')
def test_view():
    return template('upData')
    
# Tuncation tab, this tab is for tuncating the data
@get('/doTrunc')
def doTrunc():
    return template('doTrunc')

# Extraction tab, first this tab posts to /getReflection and then displays the extracted
# refractive index and the timedomain data
@get('/doExtraction')
def doExtraction():
    return template('doExtraction')

# Form post to upload the time domain data that then passes it to the extracion code
# this either returns the extracted refractive index or null if it fails
@post('/getRefractive')
def getRefractive():
    xRef = request.forms.get('xRVal')
    yRef = request.forms.get('yRVal')
    xSam = request.forms.get('xSVal')
    ySam = request.forms.get('ySVal')
    thickness = request.forms.get('thick')

    xRef = map(float, xRef.split(',')) 
    yRef = map(float, yRef.split(',')) 
    xSam = map(float, xSam.split(',')) 
    ySam = map(float, ySam.split(',')) 
    
    thickness = float(thickness)
    
    try:
        (xRVal2, yRVal2, xSVal2, ySVal2, nReal, nImag, f, start, end) = extraction.getRefrac(xRef, yRef, xSam, ySam, thickness, 1)

        f = ','.join(str(x) for x in f)    
        nI = ','.join(str(x) for x in nImag)
        nR = ','.join(str(x) for x in nReal)
    
        return { 'f': f, 'nR': nR, 'nI': nI, 'start': str(start), 'end': str(end)}
    
    except:
        return

run(host='localhost', port=8081)