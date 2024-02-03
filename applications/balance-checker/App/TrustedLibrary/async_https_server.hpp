#pragma once

#include <boost/asio.hpp>
#include <boost/asio/ssl.hpp>
#include <boost/beast.hpp>
#include <boost/beast/core/detail/base64.hpp>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
namespace asio = boost::asio;
namespace ssl = asio::ssl;
namespace beast = boost::beast;
namespace http = beast::http;
using tcp = asio::ip::tcp;

class HttpSession : public std::enable_shared_from_this<HttpSession> {
  asio::ssl::stream<beast::tcp_stream> stream_;
  beast::flat_buffer buffer_;
  http::request<http::string_body> request_;
  std::function<void(http::request<http::string_body>&,
                     http::response<http::string_body>&)>
      request_handler_;

 public:
  HttpSession(tcp::socket socket, ssl::context& ctx,
              std::function<void(http::request<http::string_body>&,
                                 http::response<http::string_body>&)>
                  handler)
      : stream_(std::move(socket), ctx), request_handler_(handler) {}

  // Start the asynchronous operation
  void start() { doHandshake(); }

 private:
  void doHandshake() {
    auto self = shared_from_this();
    stream_.async_handshake(ssl::stream_base::server,
                            [self](beast::error_code ec) {
                              if (!ec) {
                                self->readRequest();
                              }
                            });
  }

  void readRequest() {
    auto self = shared_from_this();
    http::async_read(stream_, buffer_, request_,
                     [self](beast::error_code ec, std::size_t) {
                       if (!ec) {
                         self->handleRequest();
                       }
                     });
  }

  void handleRequest() {
    http::response<http::string_body> response{http::status::ok,
                                               request_.version()};
    response.set(http::field::server, "Boost.Beast HTTPS");
    if (request_handler_) {
      request_handler_(request_, response);
    }
    response.prepare_payload();

    auto self = shared_from_this();
    http::async_write(
        stream_, response, [self](beast::error_code ec, std::size_t) {
          self->stream_.async_shutdown([](beast::error_code ec) {});
        });
  }
};

class HttpsServer {
  asio::io_context& ioc_;
  tcp::acceptor acceptor_;
  ssl::context& ssl_ctx_;
  std::function<void(const http::request<http::string_body>&,
                     http::response<http::string_body>&)>
      request_handler_;

 public:
  HttpsServer(asio::io_context& ioc, tcp::endpoint endpoint,
              ssl::context& ssl_ctx,
              std::function<void(const http::request<http::string_body>&,
                                 http::response<http::string_body>&)>
                  handler)
      : ioc_(ioc),
        acceptor_(asio::make_strand(ioc)),
        ssl_ctx_(ssl_ctx),
        request_handler_(std::move(handler)) {
    beast::error_code ec;

    acceptor_.open(endpoint.protocol(), ec);
    if (ec) {
      throw beast::system_error{ec};
    }

    acceptor_.set_option(asio::socket_base::reuse_address(true), ec);
    if (ec) {
      throw beast::system_error{ec};
    }

    acceptor_.bind(endpoint, ec);
    if (ec) {
      throw beast::system_error{ec};
    }

    acceptor_.listen(asio::socket_base::max_listen_connections, ec);
    if (ec) {
      throw beast::system_error{ec};
    }

    accept();
  }

 private:
  void accept() {
    // Accept connections and create HttpSession objects...
  }

  void onAccept(beast::error_code ec, tcp::socket socket) {
    if (!ec) {
      std::make_shared<HttpSession>(std::move(socket), ssl_ctx_,
                                    request_handler_)
          ->start();
    }

    accept();
  }
  // Other member functions...
};

int startServer() {
  try {
    // IO context
    asio::io_context ioc{1};

    // SSL context
    asio::ssl::context sslContext{asio::ssl::context::tlsv12};
    sslContext.set_options(
        asio::ssl::context::default_workarounds | asio::ssl::context::no_sslv2 |
        asio::ssl::context::no_sslv3 | asio::ssl::context::single_dh_use);

    // Load your server certificate and private key
    sslContext.use_certificate_chain_file("server.crt");
    sslContext.use_private_key_file("server.key",
                                    boost::asio::ssl::context::pem);

    // Define the endpoint (e.g., localhost on port 8080)
    auto const address = asio::ip::make_address("127.0.0.1");
    unsigned short port = 8080;

    // Instantiate and run the server
    HttpsServer server(ioc, {address, port}, sslContext,
                       [](auto& req, auto& res) {
                         res.body() = "Hello from the HTTPS server!";
                         res.prepare_payload();
                       });

    // Run the IO context to start the server
    ioc.run();
  } catch (std::exception const& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}
