#! /usr/bin/python3
import urllib, requests, lxml.html, re
#from bs4 import BeautifulSoup4

url = "http://www.cazy.org/eA.html" #http://www.cazy.org/eA.html"

htmltext = urllib.request.urlopen(url).read()
page = lxml.html.parse(htmltext)
#links = re.findall( b'href="\S*"', htmltext)
print(page)

'''
class Crawlur(object):
	"""docstring for Crawlur"""
	def __init__(self, arg):
		self.arg = arg
'''	