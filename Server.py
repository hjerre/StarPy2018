from http.server import HTTPServer
import StarPyHandler

server_address = ('127.0.0.1', 8080)

httpd = HTTPServer( server_address, StarPyHandler.StarPyHandler )
print('running server...')
httpd.serve_forever()