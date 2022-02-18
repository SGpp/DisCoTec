#include "NetworkUtils.hpp"

namespace combigrid{

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
  struct sockaddr_in servAddr;
  struct hostent *server;
  sockfd_ = socket(AF_INET, SOCK_STREAM, 0);
  if (sockfd_ < 0) {
    perror("ClientSocket::init() opening client socket failed");
    return false;
  }
  server = gethostbyname(remoteHost_.c_str());
  if (server == NULL) {
    perror(("ClientSocket::init() no such host " + remoteHost_).c_str());
    return false;
  }
  servAddr.sin_family = AF_INET;
  bcopy((char*) server->h_addr, (char*)& servAddr.sin_addr.s_addr,
      static_cast<size_t>(server->h_length));
  servAddr.sin_port = static_cast<uint16_t>(htons(remotePort_));
  int connStat = connect(sockfd_, (struct sockaddr*) &servAddr,
      sizeof(servAddr));
  if (connStat < 0) {
    perror("ClientSocket::init() connect failed");
    return false;
  }
  this->initialized_ = true;
  return true;
}

static int getSockType(int sockfd) {
  int socktype;
  socklen_t optlen = sizeof(socktype);
  getsockopt(sockfd, SOL_SOCKET, SO_TYPE, &socktype, &optlen);
  return socktype;
}

bool ClientSocket::sendall(const std::string& mesg) const {
  assert(isInitialized() && "Client Socket not initialized");
  assert(mesg.size() > 0);
  ssize_t sent = -1;
  size_t total = 0;
  const size_t& len = mesg.size();
  while (total < len) {
    sent = send(sockfd_, mesg.data() + total, len - total, 0);
    if (sent <= 0)
      break;
    total += static_cast<size_t>(sent);
  }
  switch (sent) {
    case 0:
      std::cerr << "ClientSocket::sendall() failed receiver terminated too early" << std::endl;
      return false;
    case -1:
      perror("ClientSocket::sendall() failed");
      return false;
    default:
      return true;
  }
}

bool ClientSocket::sendall(const char* buff, size_t len) const {
  assert(isInitialized() && "Client Socket not initialized");
  assert(len > 0);
  ssize_t sent = -1;
  size_t total = 0;
  while (total < len) {
    sent = send(sockfd_, &buff[total], len - total, 0);
    if (sent <= 0)
      break;
    total += static_cast<size_t>(sent);
  }
  switch (sent) {
    case 0:
      std::cerr << "ClientSocket::sendall() failed receiver terminated too early" << std::endl;
      return false;
    case -1:
      perror("ClientSocket::sendall() failed");
      return false;
    default:
      return true;
  }
}


bool ClientSocket::sendallPrefixed(const std::string& mesg) const {
  assert(isInitialized() && "Client Socket not initialized");
  assert(mesg.size() > 0);
  std::string lenstr = std::to_string(mesg.size()) + "#";
  return this->sendall(lenstr + mesg);
}

bool ClientSocket::sendallPrefixed(const char* const buff, size_t len) const {
  assert(isInitialized() && "Client Socket not initialized");
  assert(len > 0);
  std::string lenstr = std::to_string(len) + "#";
  bool ok = this->sendall(lenstr);
  if (not ok)
    return false;
  return this->sendall(buff, len);
}

bool ClientSocket::recvall(std::string& mesg, size_t len, int flags) const {
  assert(isInitialized() && "Client Socket not initialized");
  assert(len > 0);
  //std::cout << "Trying to receive " << len << "Bytes" << std::endl;
  ssize_t recvd = -1;
  std::unique_ptr<char[]> buff(new char[len]);
  auto ret = this->recvall(buff.get(), len, flags);
  mesg = std::string(buff.get(), len);
  return ret;
}

bool ClientSocket::recvallBinaryToFile(const std::string& filename, size_t len,
    size_t chunksize, int flags) const
{
  assert(isInitialized() && "Client Socket not initialized");
  assert(len > 0);
  // receive endianness first
  std::ofstream file(filename, std::ofstream::binary);
  bool success = false;

  char temp;
  success = recvall(&temp, 1);
  if (!success) {
    file.close();
    return false;
  }
  bool endianness = temp;

  ssize_t recvd = -1;
  size_t total = 0;
  std::unique_ptr<char[]> buff(new char[chunksize]);
  while (total < len) {
    recvd = recv(sockfd_, buff.get(), chunksize, 0);
    if (recvd <= 0)
      break;
    total += static_cast<size_t>(recvd);
    file.write(buff.get(), static_cast<std::streamsize>(recvd));
    file.seekp(static_cast<std::streamoff>(total)); // append to file next round
  }
  switch (recvd) {
    case 0:
      std::cerr << "ClientSocket::recvallBinaryToFile() failed, sender terminated too early" << std::endl;
      file.close();
      return false;
    case -1:
      perror("ClientSocket::recvallBinaryToFile() failed");
      file.close();
      return false;
    default:
      if (endianness != NetworkUtils::isLittleEndian()) {
        // TODO correct endiannes in file
      }
      file.close();
      return true;
  }
}

bool ClientSocket::recvall(char* buff, size_t len, int flags) const {
  assert(isInitialized() && "Client Socket not initialized");
  assert(len > 0);
  ssize_t recvd = -1;
  size_t total = 0;
  while (total < len) {
    recvd = recv(sockfd_, &buff[total], len-total, flags);
    if (recvd <= 0)
      break;
    total += static_cast<size_t>(recvd);
  }
  switch (recvd) {
    case 0:
      std::cerr << "ClientSocket::recvall() failed, sender terminated too early" << std::endl;
      return false;
    case -1:
      perror("ClientSocket::recvall() failed");
      return false;
    default:
      return true;
  }
}

bool ClientSocket::recvLength(size_t & length, int flags) const {
  // receive length string, suffixed with '#'
  ssize_t n = -1;
  std::string lenstr = "";
  char temp = ' ';
  do {
    n = recv(sockfd_, &temp, 1, flags);
    if (n <= 0)
      break;
    lenstr += temp;
  } while (temp != '#');

  switch (n) {
    case 0:
      std::cerr << "ClientSocket::recvLength() failed, sender terminated too early" << std::endl;
      return false;
    case -1:
      perror("ClientSocket::recvLength() failed");
      return false;
    default:
      lenstr.pop_back();
      assert(NetworkUtils::isInteger(lenstr) && "Received length is not a number");
      length = (size_t) std::stoi(lenstr);
      assert(length > 0);
      return true;
  }
}

bool ClientSocket::recvallPrefixed(std::unique_ptr<char[]>& buff, size_t& len, int flags) const {
  assert(isInitialized() && "Client Socket not initialized");
  if(recvLength(len, flags)) {
      // receive data
      buff.reset(new char[len]);
      return recvall(buff.get(), len);
  }
  return false;
}

bool ClientSocket::recvallPrefixed(std::string& mesg, int flags) const {
  assert(isInitialized() && "Client Socket not initialized");
  mesg.clear();
  size_t len;
  if(recvLength(len, flags)) {
      return recvall(mesg, len);
  }
  return false;
}

bool ClientSocket::isReadable(int timeoutSec) const {
  assert(isInitialized() && "Client Socket not initialized");

  struct timeval tv;
  fd_set readfds;

  FD_ZERO(&readfds);
  FD_SET(sockfd_, &readfds);

  if (timeoutSec > -1) {
    tv.tv_sec = timeoutSec;
    tv.tv_usec = 0;
    select(sockfd_+1, &readfds, NULL, NULL, &tv);
  }
  else {
    select(sockfd_+1, &readfds, NULL, NULL, NULL);
  }
  return FD_ISSET(sockfd_, &readfds);
}

bool ClientSocket::isWriteable(int timeoutSec) const {
  assert(isInitialized() && "Client Socket not initialized");

  struct timeval tv;
  fd_set writefds;

  tv.tv_sec = timeoutSec;
  tv.tv_usec = 0;

  FD_ZERO(&writefds);
  FD_SET(sockfd_, &writefds);

  if (timeoutSec > -1)
    select(sockfd_+1, NULL, &writefds, NULL, &tv);
  else
    select(sockfd_+1, NULL, &writefds, NULL, NULL);
  return FD_ISSET(sockfd_, &writefds);
}

std::string ClientSocket::getRemoteHost() const {
  return remoteHost_;
}

int ClientSocket::getRemotePort() const {
  return remotePort_;
}

sockaddr_in ClientSocket::getSockaddrIn() const {
  struct sockaddr_in addr;
  socklen_t socklen = sizeof(addr);
  bzero(&addr, socklen);
  addr.sin_family = AF_INET;
  addr.sin_port = htons(remotePort_);
  if (inet_aton(remoteHost_.c_str(), &addr.sin_addr) <= 0) {
    throw std::runtime_error("not a valid IPaddress " + remoteHost_);
  }
  return addr;
}

ServerSocket::ServerSocket() : Socket(), port_(0) {
}

ServerSocket::ServerSocket(const unsigned short port) : Socket(), port_(port) {
}

ServerSocket::~ServerSocket() {
  close(sockfd_);
}

bool ServerSocket::init() {
  struct sockaddr_in servAddr;

  sockfd_ = socket( AF_INET, SOCK_STREAM, 0 );
  if (sockfd_ < 0) {
    perror("ServerSocket::init() opening server socket failed");
    return false;
  }

  bzero((char*) &servAddr, sizeof(servAddr));
  servAddr.sin_family = AF_INET;
  servAddr.sin_port = static_cast<uint16_t>(htons(port_));
  servAddr.sin_addr.s_addr = INADDR_ANY;
  int bindstat = bind(sockfd_, (struct sockaddr*) &servAddr, sizeof(servAddr));
  if (bindstat < 0) {
    perror(("ServerSocket::init() binding to port " + std::to_string(port_) + " failed").c_str());
    return false;
  }

  int listenstat = listen(sockfd_, 1);
  if (listenstat < 0) {
    perror("ServerSocket::init() listen failed");
    return false;
  }

  // set port if determined by os
  if (this->port_ == 0) {
    socklen_t len = sizeof(servAddr);
    int stat =  getsockname(sockfd_, (struct sockaddr *)&servAddr, &len);
    if (stat < 0) {
      perror("ServerSocket::init() querying port failed");
      return false;
    } else {
      port_ = ntohs(servAddr.sin_port);
      std::cout << "Server runs on port: " << port_ << std::endl;
    }
  }
  initialized_ = true;
  return true;
}

std::shared_ptr<ClientSocket> ServerSocket::acceptClient() const {
  assert(isInitialized() && "Server Socket not initialized");
  std::shared_ptr<ClientSocket> client;

  // blocks until a new client connects
  struct sockaddr_in cliAddr;
  socklen_t cliLen = sizeof(cliAddr);
  int clientfd = accept(sockfd_, (struct sockaddr*) &cliAddr, &cliLen);
  if (clientfd < 0) {
    perror("ServerSocket::acceptClient() accept failed");
    return client;
  }
  // initialize ClientSocket
  client.reset(new ClientSocket);
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

Socket::~Socket(){
  close(sockfd_);
}

int Socket::getFileDescriptor() const {
  return this->sockfd_;
}

bool Socket::isInitialized() const {
  return initialized_ && sockfd_ > 0;
}

bool NetworkUtils::forward(const ClientSocket& sender,
    const ClientSocket& receiver,  size_t chunksize, size_t size)
{
  assert(sender.isInitialized() && "Initialize sender first");
  assert(receiver.isInitialized() && "Initialize receiver first");
  size_t totalRecvd = 0;
  ssize_t recvd = 0;
  bool sendSuccess = false;
  int sendFd = sender.getFileDescriptor();
  std::unique_ptr<char[]> buff(new char[chunksize]);
  while (totalRecvd < size)
  {
#ifdef DEBUG_OUTPUT
    std::cout << "." << std::flush;
#endif
    size_t remaining = size - totalRecvd;
    if (remaining < chunksize)
      recvd = recv(sendFd, buff.get(), remaining, 0);
    else
      recvd = recv(sendFd, buff.get(), chunksize, 0);
    switch (recvd) {
      case 0:
        std::cerr << "NetworkUtils::forward() sender terminated too early" << std::endl;
        return false;
      case -1:
        perror("NetworkUtils::forward() unexpected fail of sender");
        return false;
      default:
        totalRecvd += static_cast<size_t>(recvd);
        // send received bytes to receiver
        sendSuccess = receiver.sendall(buff.get(), static_cast<size_t>(recvd));
        if (!sendSuccess) {
          std::cerr << "NetworkUtils::forward() unexpected fail of receiver";
          return false;
        }
    }
  }
#ifdef DEBUG_OUTPUT
  std::cout << std::endl;
#endif
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

/*
 * Splits a string into tokens along the delimiter c
 */
void NetworkUtils::split(const std::string& s, const char c,
   std::vector<std::string>& tokens)
{
  size_t next;
  size_t prev = 0;
  tokens.clear();
  while((next = s.find(c, prev)) != std::string::npos) {
    tokens.push_back(s.substr(prev, next-prev));
    prev = next + 1;
  }
  if (prev < s.size())
    tokens.push_back(s.substr(prev));
}

} // namespace combigrid