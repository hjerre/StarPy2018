from http.server import BaseHTTPRequestHandler
from urllib.parse import urlparse


class StarPyHandler( BaseHTTPRequestHandler ):

    def do_GET(self):
        self.do_POST()

    def do_POST(self):
        print("I Handler med do_Post")
        elements = urlparse(self.path)
        sender = self.client_address
        request = elements.path
        print("Query   : ", elements.query)
        func = elements.query.split("&")
        for elem in func:
            tmp = elem.split("=")
            print(tmp[0], " : ", tmp[1])
