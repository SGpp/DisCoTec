#include "NetworkUtils.hpp"

ClientSocket::ClientSocket(const std::string& host, const int port)
  : Socket(), remotePort_(port), remoteHost_(host){
}

ClientSocket::ClientSocket() : Socket() {
  initialized_ = false;
}

ClientSocket::~ClientSocket() {
  close(sockfd_);
}

bool ClientSocket::init() {
  assert(!isInitialized() && "Client is already initialized");
  int err;
  struct sockaddr_in servAddr;
  struct hostent *server;
  sockfd_ = socket(AF_INET, SOCK_STREAM, 0);
  err = errno;
  if (sockfd_ < 0) {
    std::cerr << "Opening client socket failed: " << err << std::endl;
    return false;
  }
  server = gethostbyname(remoteHost_.c_str());
  err = errno;
  if (server == NULL) {
    std::cerr << "No such host " << remoteHost_ << ": " << err
      << std::endl;
    return false;
  }
  servAddr.sin_family = AF_INET;
  bcopy((char*) server->h_addr, (char*)& servAddr.sin_addr.s_addr,
      server->h_length);
  servAddr.sin_port = htons(remotePort_);
  int connStat = connect(sockfd_, (struct sockaddr*) &servAddr,
      sizeof(servAddr));
  err = errno;
  if (connStat < 0) {
    std::cerr << "Connect failed: " << err << std::endl;
    return false;
  }
  this->initialized_ = true;
  return true;
}

bool ClientSocket::sendall(const std::string& mesg) const {
  assert(isInitialized() && "Client Socket not initialized");
  int err;
  long n;
  size_t total = 0;
  const size_t& len = mesg.size();
  while (total < len) {
    n = send(sockfd_, mesg.data() + total, len - total, 0);
    err = errno;
    if (n == -1)
      break;
    total += n;
  }
  if ( n == -1) {
    std::cerr << "Send failed: " << std::to_string(err) << std::endl;
    return false;
  }
  return true;
}

bool ClientSocket::sendall(const char* buf, size_t len) const {
  assert(isInitialized() && "Client Socket not initialized");
  int err;
  long n;
  size_t total = 0;
  while (total < len) {
    n = send(sockfd_, &buf[total], len - total, 0);
    err = errno;
    if (n == -1)
      break;
    total += n;
  }
  if ( n == -1) {
    std::cerr << "Send failed: " << std::to_string(err) << std::endl;
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
    n = recv(sockfd_, &buff[total], len-total, flags);
    err = errno;
    if (n == -1)
      break;
    total += n;
  }
  if ( n == -1) {
    std::cerr << "Receive failed: " << std::to_string(err) << std::endl;
    return false;
  }
  mesg = std::string(buff, len);
  delete[] buff;
  return true;
}

bool ClientSocket::recvallBinaryToFile(const std::string& filename, size_t len, size_t chunksize) const {
  // receive endianness first
  assert(isInitialized() && "Client Socket not initialized");
  std::ofstream file(filename, std::ofstream::binary);
  bool success = false;
  //recv endianness first
  char* endianness = nullptr;
  success = recvall(endianness, 1);
  if (!success) {
    file.close();
    return false;
  }

  int err;
  long n;
  size_t total = 0;
  char* buff = new char[chunksize];
  while (total < len) {
    n = recv(sockfd_, buff, chunksize, 0);
    err = errno;
    if (n == -1)
      break;
    total += static_cast<size_t>(n);
    file.write(buff, static_cast<std::streamsize>(n));
    file.seekp(static_cast<std::streamoff>(total)); // start writing at end of file
  }
  if ( n == -1) {
    std::cerr << "Receive failed: " << std::to_string(err) << std::endl;
    file.close();
    return false;
  }

  if ((bool) *endianness != NetworkUtils::isLittleEndian()) {
    // TODO correct endiannes of file
  }

  file.close();
  delete[] endianness;
  delete[] buff;
  return true;
}

bool ClientSocket::recvall(char* &buf, size_t len, int flags) const {
  assert(isInitialized() && "Client Socket not initialized");
  int err;
  long n;
  buf = new char[len];
  size_t total = 0;
  while (total < len) {
    n = recv(sockfd_, &buf[total], len-total, flags);
    err = errno;
    if (n == -1)
      break;
    total += n;
  }
  if ( n == -1) {
    std::cerr << "Receive failed: " << std::to_string(err) << std::endl;
    return false;
  }
  return true;
}

bool ClientSocket::recvallPrefixed(std::string& mesg, int flags) const {
  assert(isInitialized() && "Client Socket not initialized");
  int err;
  long n;
  std::string lenstr = "";
  char temp = ' ';
  do {
    n = recv(sockfd_, &temp, 1, flags);
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

  FD_ZERO(&readfds);
  FD_SET(sockfd_, &readfds);

  if (timeout > -1)
    select(sockfd_+1, &readfds, NULL, NULL, &tv);
  else
    select(sockfd_+1, &readfds, NULL, NULL, NULL);
  return FD_ISSET(sockfd_, &readfds);
}

bool ClientSocket::isWriteable(int timeout) const {
  struct timeval tv;
  fd_set writefds;

  tv.tv_sec = timeout;
  tv.tv_usec = 0;

  FD_ZERO(&writefds);
  FD_SET(sockfd_, &writefds);

  if (timeout > -1)
    select(sockfd_+1, NULL, &writefds, NULL, &tv);
  else
    select(sockfd_+1, NULL, &writefds, NULL, NULL);
  return FD_ISSET(sockfd_, &writefds);
}

ServerSocket::ServerSocket(const unsigned short port) : Socket(), port_(port) {
  init();
}

ServerSocket::ServerSocket() : Socket(), port_(0) {
  init();
}

ServerSocket::~ServerSocket() {
  close(sockfd_);
}

bool ServerSocket::init() {
  int err;
  struct sockaddr_in servAddr;

  sockfd_ = socket( AF_INET, SOCK_STREAM, 0 );
  err = errno;
  if (!isInitialized()) {
    std::cerr << "Opening server socket failed: " << err << std::endl;
    return false;
  }

  bzero((char*) &servAddr, sizeof(servAddr));
  servAddr.sin_family = AF_INET;
  servAddr.sin_port = htons(port_);
  servAddr.sin_addr.s_addr = INADDR_ANY;
  int bindstat = bind(sockfd_, (struct sockaddr*) &servAddr, sizeof(servAddr));
  err = errno;
  if (bindstat < 0) {
    std::cerr << "Binding to port " << port_ << " failed: " << err
      << std::endl;
    return false;
  }

  int listenstat = listen(sockfd_, 1);
  err = errno;
  if (listenstat < 0) {
    std::cerr << "Listen failed: " << err << std::endl;
    return false;
  }

  // set port if determined by os (port == 0)
  if (this->port_ == 0) {
    socklen_t len = sizeof(servAddr);
    int stat =  getsockname(sockfd_, (struct sockaddr *)&servAddr, &len);
    err = errno;
    if (stat < 0) {
      std::cerr << "Querying port failed: " << err
        << std::endl;
      return false;
    } else {
      port_ = ntohs(servAddr.sin_port);
    }
  }
  initialized_ = true;
  return true;
}

ClientSocket* ServerSocket::acceptClient() const {
  assert(isInitialized() && "Server Socket not initialized");
  ClientSocket* client = new ClientSocket();
  int err;
  // blocks until a new client connects
  struct sockaddr_in cliAddr;
  socklen_t cliLen = sizeof(cliAddr);
  int clientfd = accept(sockfd_, (struct sockaddr*) &cliAddr, &cliLen);
  err = errno;
  if (clientfd < 0) {
    std::cerr << "Accept failed: " << err << std::endl;
    return client;
  }
  // initialize ClientSocket
  unsigned short port = cliAddr.sin_port;
  std::string host = std::string(inet_ntoa(cliAddr.sin_addr));
  client->remoteHost_ = host;
  client->remotePort_ = port;
  client->sockfd_ = clientfd;
  client->initialized_ = true;
  return client;
}

int ServerSocket::getPort() {
  return this->port_;
}

Socket::Socket() : sockfd_(-1), initialized_(false) {
}

int Socket::getFileDescriptor() const {
  return this->sockfd_;
}

bool Socket::isInitialized() const {
  return initialized_ && sockfd_ > 0;
}

bool NetworkUtils::tunnel(const ClientSocket* sender, const ClientSocket* receiver,
    unsigned int buffsize) {
  assert(sender->isInitialized() && "Initialize sender first");
  assert(receiver->isInitialized() && "Initialize receiver first");
  int err;
  int totalRecvd = 1;
  int sendFd = sender->getFileDescriptor();
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
    bool success = receiver->sendall(std::string(buff, totalRecvd));
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

/*
 * Checks if callers system is little-endian byte order.
 */
bool NetworkUtils::isLittleEndian() {
  int x = 1;
  return (*(char*)&x == 1);
}

template<typename T>
  bool sendBinary(const std::vector<T>& buf, const ClientSocket& socket)
  {
    if (!socket.isInitialized())
      return false;
    bool success = false;
    size_t size = buf.size() * sizeof(T);
    char endianness = NetworkUtils::isLittleEndian();
    // send endianness header
    success = socket.sendall(&endianness, 1);
    if (!success)
      return false;
    success = socket.sendall((char*)buf.data(), size);
    if (!success)
      return false;
    return true;
  }

// TODO Optimize memory usage (hopefully this wont cause thrashing and os recognizes the linear access)
template<typename T>
  bool recvBinary(std::vector<T>& buf, size_t len, const ClientSocket& socket)
  {
    if (!socket.isInitialized())
      return false;
    bool success = false;
    //recv endianness first
    char* endianess = nullptr;
    success = socket.recvall(endianess, 1);
    if (!success)
      return false;
    //recv data
    char* raw = nullptr;
    size_t rawSize = len*sizeof(T);
    success = socket.recvall(raw, rawSize);
    if (!success)
      return false;
    buf.clear();
    buf.reserve(len);
    if ((bool) *endianess == NetworkUtils::isLittleEndian()) { //same endianness
      T* values = reinterpret_cast<T*>(raw);
      // copy values to buf
      for (size_t i = 0; i < len; i++) { // linear access to poss. huge arrays
        buf.push_back(values[i]);
      }
    } else { // different endianness
      char* rightOrder = new char[sizeof(T)];
      // extract values from raw buffer in right order
      // reverts byte order for each sizeof(T) big block
      for (size_t i = 0; i < len; i+=sizeof(T)) { // somehow linear access
        for (size_t j = sizeof(T)-1; j >= 0; j--) {
          rightOrder[sizeof(T)-1 - j] = raw[i+j];
          buf.push_back( reinterpret_cast<T &>(rightOrder) );
        }
      }
      delete[] rightOrder;
    }
    delete[] endianess;
    delete[] raw;
    return true;
  }
