#include "NetworkUtils.hpp"

ClientSocket::ClientSocket(const std::string& host, const int port)
  : Socket(), remotePort(port), remoteHost(host){
}

bool ClientSocket::init() {
  assert(!isInitialized() && "Client is already initialized");
  int err;
  struct sockaddr_in servAddr;
  struct hostent *server;
  sockfd = socket(AF_INET, SOCK_STREAM, 0);
  err = errno;
  if (sockfd < 0) {
    std::cerr << "Opening client socket failed: " << err << std::endl;
    return false;
  }
  server = gethostbyname(this->remoteHost.c_str());
  err = errno;
  if (server == NULL) {
    std::cerr << "No such host " << this->remoteHost << ": " << err
      << std::endl;
    return false;
  }
  servAddr.sin_family = AF_INET;
  bcopy((char*) server->h_addr, (char*)& servAddr.sin_addr.s_addr,
      server->h_length);
  servAddr.sin_port = htons(this->remotePort);
  int connStat = connect(sockfd, (struct sockaddr*) &servAddr,
      sizeof(servAddr));
  err = errno;
  if (connStat < 0) {
    std::cerr << "Connect failed: " << err << std::endl;
    return false;
  }
  this->initialized = true;
  return true;
}

bool ClientSocket::sendall(const std::string& mesg) const {
  assert(isInitialized() && "Client Socket not initialized");
  int err;
  long n;
  unsigned int total = 0;
  const size_t& len = mesg.size();
  while (total < len) {
    n = send(sockfd, mesg.data() + total, len - total, 0);
    err = errno;
    if (n == -1)
      break;
    total += n;
  }
  if ( n == -1) {
    std::cerr << "Send failed: " << std::to_string(err) <<
      std::endl;
    return false;
  }
  return true;
}

bool ClientSocket::sendallPrefixed(const std::string& mesg) const {
  assert(isInitialized() && "Client Socket not initialized");
  std::string lenstr = std::to_string(mesg.size()) + "#";
  return this->sendall(lenstr + mesg);
}

bool ClientSocket::recvall(std::string& mesg, size_t len, int flags) const {
  assert(isInitialized() && "Client Socket not initialized");
  int err;
  long n;
  char* buff = new char[len];
  size_t total = 0;
  while (total < len) {
    n = recv(sockfd, &buff[total], len-total, flags);
    err = errno;
    if (n == -1)
      break;
    total += n;
  }
  if ( n == -1) {
    std::cerr << "Receive failed: " << std::to_string(err);
    return false;
  }
  mesg = std::string(buff, len);
  delete[] buff;
  return true;
}

bool ClientSocket::recvall(std::stringstream& mesg, size_t len, int flags) const {
  assert(isInitialized() && "Client Socket not initialized");
  int err;
  long n;
  char* buff = new char[len];
  size_t total = 0;
  while (total < len) {
    n = recv(sockfd, &buff[total], len-total, flags);
    err = errno;
    if (n == -1)
      break;
    total += n;
  }
  if ( n == -1) {
    std::cerr << "Receive failed: " << std::to_string(err);
    return false;
  }
  mesg = std::string(buff, len);
  return true;
}

bool ClientSocket::recvallPrefixed(std::string& mesg, int flags) const {
  assert(isInitialized() && "Client Socket not initialized");
  int err;
  long n;
  std::string lenstr = "";
  char temp = ' ';
  do {
    n = recv(sockfd, &temp, 1, flags);
    err = errno;
    if (n == -1)
      break;
    lenstr += temp;
  } while (temp != '#');
  if (n == -1) {
    std::cerr << "Receive of length failed: " << std::to_string(err);
    return false;
  }
  lenstr.pop_back();
  assert(NetworkUtils::isInteger(lenstr) && "Received length is not a number");
  size_t len = (size_t) std::stoi(lenstr);
  return recvall(mesg, len);
}


bool ClientSocket::isReadable(int timeout) const {
  struct timeval tv;
  fd_set readfds;

  tv.tv_sec = timeout;
  tv.tv_usec = 0;

  FD_ZERO(&readfds); // clear readfds
  FD_SET(sockfd, &readfds);

  if (timeout > -1)
    select(sockfd+1, &readfds, NULL, NULL, &tv);
  else
    select(sockfd+1, &readfds, NULL, NULL, NULL);
  return FD_ISSET(sockfd, &readfds);
}

bool ClientSocket::isWriteable(int timeout) const {
  struct timeval tv;
  fd_set writefds;

  tv.tv_sec = timeout;
  tv.tv_usec = 0;

  FD_ZERO(&writefds); // clear writefds
  FD_SET(sockfd, &writefds);

  if (timeout > -1)
    select(sockfd+1, NULL, &writefds, NULL, &tv);
  else
    select(sockfd+1, NULL, &writefds, NULL, NULL);
  return FD_ISSET(sockfd, &writefds);
}

ServerSocket::ServerSocket(const unsigned short port) : Socket(), port(port) {
    init();
}

bool ServerSocket::init() {
  int err;
  struct sockaddr_in servAddr;
  sockfd = socket( AF_INET, SOCK_STREAM, 0 );
  err = errno;
  if (!isInitialized()) {
    std::cerr << "Opening server socket failed: " << err << std::endl;
    return false;
  }
  bzero((char*) &servAddr, sizeof(servAddr));
  servAddr.sin_family = AF_INET;
  servAddr.sin_port = htons(this->port);
  servAddr.sin_addr.s_addr = INADDR_ANY;
  int bindstat = bind(sockfd, (struct sockaddr*) &servAddr, sizeof(servAddr));
  err = errno;
  if (bindstat < 0) {
    std::cerr << "Binding to port " << this->port << " failed: " << err
      << std::endl;
    return false;
  }
  this->initialized = true;
  return true;
}

ClientSocket ServerSocket::acceptClient() const {
  assert(isInitialized() && "Server Socket not initialized");
  ClientSocket client("", 0);
  int err;
  int listenstat = listen(sockfd, 1);
  err = errno;
  if (listenstat < 0) {
    std::cerr << "Listen failed: " << err << std::endl;
    return client;
  }
  // blocks until a new client connects
  struct sockaddr_in cliAddr;
  socklen_t cliLen = sizeof(cliAddr);
  int clientfd = accept(this->sockfd, (struct sockaddr*) &cliAddr, &cliLen);
  err = errno;
  if (clientfd < 0) {
    std::cerr << "Accept failed: " << err << std::endl;
    return client;
  }

  // initialize ClientSocket
  unsigned short port = cliAddr.sin_port;
  // hopefully this is null byte terminated
  std::string host = std::string(inet_ntoa(cliAddr.sin_addr));
  client.remoteHost = host;
  client.remotePort = port;
  client.sockfd = clientfd;
  client.initialized = true;
  return client;
}

Socket::Socket() : sockfd(-1), initialized(false) {
}

int Socket::getFileDescriptor() const {
  return this->sockfd;
}

bool Socket::isInitialized() const {
  return initialized && sockfd > 0;
}

bool NetworkUtils::tunnel(ClientSocket sender, ClientSocket receiver,
    unsigned int buffsize) {
  assert(sender.isInitialized() && "Initialize sender first");
  assert(receiver.isInitialized() && "Initialize receiver first");
  int err;
  int totalRecvd = 1;
  int sendFd = sender.getFileDescriptor();
  char* buff = new char[buffsize];
  while (totalRecvd < buffsize) {
    totalRecvd = recv(sendFd, buff, buffsize, 0);
    err = errno;
    if (totalRecvd == -1) {
        std::cerr << "Unexpected fail of sender" << err;
        delete[] buff;
        return false;
    }
    /* can be improved by sending and receiving asynchronous in
     * separate threads*/
    bool success = receiver.sendall(std::string(buff, totalRecvd));
    if (!success) {
      std::cerr << "Unexpected fail of receiver";
      delete[] buff;
      return false;
    }
  }
  delete[] buff;
  return true;
}

/*
 * Checks if a given string represents a decimal integer.
 */
bool NetworkUtils::isInteger(const std::string& s) {
  std::string::const_iterator it = s.begin();
  if (*it == '+' || *it == '-') ++it;
  while(it != s.end() && std::isdigit(*it)) ++it;
  return !s.empty() && it == s.end();
}
