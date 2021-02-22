import http.server
import socketserver

PORT = 8080

class quietServer(http.server.SimpleHTTPRequestHandler):
    def log_message(self, format, *args):
        pass

with socketserver.TCPServer(("", PORT), quietServer) as httpd:
    httpd.serve_forever()