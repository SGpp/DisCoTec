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
  long n = -1;
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
  long n = -1;
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

template<typename FG_ELEMENT>
bool ClientSocket::sendallBinary(std::vector<FG_ELEMENT>& buf, int flags) const {
  char* rawBuf = reinterpret_cast<char*>(buf.data());
  size_t rawSize = buf.data() * sizeof(FG_ELEMENT);

  char isLittleEndian = static_cast<char>(NetworkUtils::isLittleEndian());
  bool success = sendall(&isLittleEndian, 1);
  if (!success) 
    return false;
  return sendall(rawBuf, rawSize);
}

// Sends the given buffer in binary form. The first byte sent specifies
// if the system is little-endian. The next bytes specify the length of
// the following buffer and are sperated by an # from the following raw data. 
template<typename FG_ELEMENT>
bool ClientSocket::sendallBinaryPrefixed(std::vector<FG_ELEMENT>& buf, int flags) const {
  char* rawBuf = reinterpret_cast<char*>(buf.data());
  size_t rawSize = buf.data() * sizeof(FG_ELEMENT);

  char isLittleEndian = static_cast<char>(NetworkUtils::isLittleEndian());
  bool success = sendall(&isLittleEndian, 1);
  if (!success)
    return false;
  return sendallPrefixed(rawBuf, rawSize);
}

bool ClientSocket::sendallPrefixed(const std::string& mesg) const {
  assert(isInitialized() && "Client Socket not initialized");
  std::string lenstr = std::to_string(mesg.size()) + "#";
  return this->sendall(lenstr + mesg);
}

bool ClientSocket::sendallPrefixed(const char* buf, size_t len) const {
  assert(isInitialized() && "Client Socket not initialized");
  std::string lenstr = std::to_string(len) + "#";
  bool ok = this->sendall(lenstr);
  if (not ok)
    return ok;
  return this->sendall(buf, len);
}

bool ClientSocket::recvall(std::string& mesg, size_t len, int flags) const {
  assert(isInitialized() && "Client Socket not initialized");
  int err;
  long n = -1;
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

bool ClientSocket::recvallBinaryToFile(const std::string& filename, size_t len,
    size_t chunksize, int flags) const
{
  // receive endianness first
  assert(isInitialized() && "Client Socket not initialized");
  std::ofstream file(filename, std::ofstream::binary);
  bool success = false;

  char temp;
  success = recvall(&temp, 1);
  if (!success) {
    file.close();
    return false;
  }
  bool endianness = temp;

  int err;
  long n = -1;
  size_t total = 0;
  char* buff = new char[chunksize];
  while (total < len) {
    n = recv(sockfd_, buff, chunksize, 0);
    err = errno;
    if (n <= 0)
      break;
    total += static_cast<size_t>(n);
    file.write(buff, static_cast<std::streamsize>(n));
    file.seekp(static_cast<std::streamoff>(total)); // append to file next round
  }
  if ( n == -1) {
    std::cerr << "Receive failed: " << std::to_string(err) << std::endl;
    file.close();
    return false;
  } else if (n == 0) {
    std::cerr << "Receive failed, sender terminated too early: " <<
      std::to_string(err) << std::endl;
    return false;
  }

  if (endianness != NetworkUtils::isLittleEndian()) {
    // TODO correct endiannes in file
  }

  file.close();
  delete[] buff;
  return true;
}

bool ClientSocket::recvall(char* buf, size_t len, int flags) const {
  assert(isInitialized() && "Client Socket not initialized");
  int err;
  long n = -1;
  size_t total = 0;
  while (total < len) {
    n = recv(sockfd_, &buf[total], len-total, flags);

    err = errno;
    if (n <= 0)
      break;
    total += n;

    if (len  == 67584)
      std::cout << "received: " << total << std::endl;

  }
  if ( n == -1) {
    std::cerr << "Receive failed: " << std::to_string(err) << std::endl;
    return false;
  } else if (n == 0) {
    std::cerr << "Receive failed, sender terminated too early: " <<
      std::to_string(err) << std::endl;
    return false;
  }
  return true;
}

bool ClientSocket::recvallPrefixed(char* buf, int flags) const {
  assert(isInitialized() && "Client Socket not initialized");

  // receive length
  int err;
  long n = -1;
  std::string lenstr = "";
  char temp = ' ';
  do {
    n = recv(sockfd_, &temp, 1, flags);
    err = errno;
    if (n <= 0)
      break;
    lenstr += temp;
  } while (temp != '#');
  if (n == -1) {
    std::cerr << "Receive of length failed: " << std::to_string(err)
      << std::endl;
    return false;
  } else if (n == 0) {
    std::cerr << "Receive of length failed, sender terminated too early: " <<
      std::to_string(err) << std::endl;
    return false;
  }
  lenstr.pop_back();
  assert(NetworkUtils::isInteger(lenstr) && "Received length is not a number");
  size_t len = (size_t) std::stoi(lenstr);

  // receive data
  return recvall(buf, len);
}

bool ClientSocket::recvallPrefixed(std::string& mesg, int flags) const {
  assert(isInitialized() && "Client Socket not initialized");

  // receive length
  mesg.clear();
  int err;
  long n = -1;
  std::string lenstr = "";
  char temp = ' ';
  do {
    n = recv(sockfd_, &temp, 1, flags);
    err = errno;
    if (n <= 0)
      break;
    lenstr += temp;
  } while (temp != '#');
  if (n == -1) {
    std::cerr << "Receive of length failed: " << std::to_string(err)
      << std::endl;
    return false;
  } else if (n == 0) {
    std::cerr << "Receive of length failed, sender terminated too early: " <<
      std::to_string(err) << std::endl;
    return false;
  }
  lenstr.pop_back();
  assert(NetworkUtils::isInteger(lenstr) && "Received length is not a number");
  size_t len = (size_t) std::stoi(lenstr);

  // receive data
  return recvall(mesg, len);
}


template<typename FG_ELEMENT>
bool ClientSocket::recvallBinaryPrefixedCorrectInPlace(std::vector<FG_ELEMENT>& buff, size_t chunksize, int flags) const {
  assert(isInitialized() && "Client Socket not initialized");
  // receive endianness
  char temp = ' ';
  assert(recvall(&temp, 1)&& "Receiving Endianess failed");
  bool endianness = temp;

  // receive length
  int err;
  long n = -1;
  std::string lenstr = "";
  do {
    n = recv(sockfd_, &temp, 1, flags);
    err = errno;
    if (n <= 0)
      break;
    lenstr += temp;
  } while (temp != '#');
  if (n == -1) {
    std::cerr << "Receive of length failed: " << std::to_string(err)
      << std::endl;
    return false;
  } else if (n == 0) {
    std::cerr << "Receive of length failed, sender terminated too early: " <<
      std::to_string(err) << std::endl;
    return false;
  }
  lenstr.pop_back();
  assert(NetworkUtils::isInteger(lenstr) && "Received length is not a number");
  size_t len = (size_t) std::stoi(lenstr);

  buff.resize(len);
  char* rawBuffer = reinterpret_cast<char*>(buff.data());

  // receive binary data and correct endianness in place
  std::vector<char> tempBuff(chunksize);
  std::vector<char> overhead(sizeof(FG_ELEMENT));
  size_t total = 0;
  size_t oldSize = 0;
  size_t newSize = 0;
  while(total < len) {
    n = recv(sockfd_, tempBuff.data(), len-total, flags);

    err = errno;
    if (n <= 0)
      break;

    // tempBuff example for sizeof(FG_ELEMENT)=4:
    //
    //        recv k-1                     recv k
    //                     old     tail               new
    //  [ oo|xxxx|...|xxxx|o ] -> [ ooo|xxxx|...|xxxx|oo ]
    //                              0                  n-1
    if (total == 0) //first run
      oldSize = n % sizeof(FG_ELEMENT);
    else
      oldSize = newSize;

    if (endianness != NetworkUtils::isLittleEndian()) {
      size_t tailSize = sizeof(FG_ELEMENT) - oldSize;
      size_t rawSize = total + tailSize;
      newSize = (tailSize - n) % sizeof(FG_ELEMENT);
      size_t blocksSize = (n - tailSize - newSize);
      if (tailSize > 0) {
        // fill remaining block
        for (size_t i = 0; i < tailSize; i++)
          rawBuffer[total-oldSize-1 + i] = tempBuff[i];
      }
      if (blocksSize > 0) {
        // append correct sized blocks
        for (size_t i = rawSize; i < rawSize+n-newSize; i+=sizeof(FG_ELEMENT)) {
          for (size_t j = sizeof(FG_ELEMENT)-1; j >= 0; j--) {
            rawBuffer[i+(sizeof(FG_ELEMENT)-1-j)] = tempBuff[i-rawSize+j];
          }
        }
      }
      if (newSize > 0) {
        rawSize += n-newSize;
        // append new overhead
        for (size_t i = rawSize; i < rawSize+newSize; i++){
          size_t j = n - (i-rawSize) -1;
          rawBuffer[i] = tempBuff[j];
        }
      }
    }

    total += n;
  }
  if ( n == -1) {
    std::cerr << "Receive failed: " << std::to_string(err) << std::endl;
    return false;
  } else if (n == 0) {
    std::cerr << "Receive failed, sender terminated too early: " <<
      std::to_string(err) << std::endl;
    return false;
  }
  return true;
}

template<typename FG_ELEMENT>
bool ClientSocket::recvallBinaryPrefixed(std::vector<FG_ELEMENT>& buff, bool& endianness, int flags) const
{
  assert(isInitialized() && "Client Socket not initialized");

  // receive endianness
  char temp = ' ';
  assert(recvall(&temp, 1)&& "Receiving Endianess failed");
  endianness = temp;

  // receive binary data
  char* rawBuf = reinterpret_cast<char*>(buff.data());
  return recvallPrefixed(rawBuf, flags);

}

template<typename FG_ELEMENT>
bool ClientSocket::recvallBinary(std::vector<FG_ELEMENT>& buff, size_t len,
    bool& endianness, int flags) const {
  // receive endianness
  char temp = ' ';
  assert(recvall(&temp, 1)&& "Receiving Endianess failed");
  endianness = temp;

  // receive binary data
  char* rawBuf = reinterpret_cast<char*>(buff.data());
  return recvall(rawBuf, len*sizeof(FG_ELEMENT), flags);
}

bool ClientSocket::isReadable(int timeout) const {
  assert(isInitialized() && "Client Socket not initialized");

  struct timeval tv;
  fd_set readfds;


  FD_ZERO(&readfds);
  FD_SET(sockfd_, &readfds);

  if (timeout > -1) {
    tv.tv_sec = timeout;
    tv.tv_usec = 0;
    select(sockfd_+1, &readfds, NULL, NULL, &tv);
  }
  else {
    select(sockfd_+1, &readfds, NULL, NULL, NULL);
  }
  return FD_ISSET(sockfd_, &readfds);
}

bool ClientSocket::isWriteable(int timeout) const {
  assert(isInitialized() && "Client Socket not initialized");

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

std::string ClientSocket::getRemoteHost() const {
  return remoteHost_;
}

int ClientSocket::getRemotePort() const {
  return remotePort_;
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
  if (sockfd_ < 0) {
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

  // set port if determined by os
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

Socket::~Socket(){
  close(sockfd_);
}

int Socket::getFileDescriptor() const {
  return this->sockfd_;
}

bool Socket::isInitialized() const {
  return initialized_ && sockfd_ > 0;
}

bool NetworkUtils::forward(const ClientSocket* sender,
    const ClientSocket* receiver, size_t size, size_t chunksize)
{
  assert(sender->isInitialized() && "Initialize sender first");
  assert(receiver->isInitialized() && "Initialize receiver first");
  int err;
  size_t totalRecvd = 0;
  long recvd = 0;
  bool sendSuccess = false;
  int sendFd = sender->getFileDescriptor();
  char* buff = new char[chunksize];
  if (size != 0) {
    std::cout << "Start tunneling of " << size << " Bytes with same Endianess:" << std::endl;
    while (totalRecvd < size)
    {
      std::cout << ".";
      // receive a max of chunksize bytes
      recvd = recv(sendFd, buff, chunksize, 0);
      err = errno;
      if (recvd == -1) {
          std::cerr << "Unexpected fail of sender" << err;
          delete[] buff;
          return false;
      } else if (recvd == 0) {
          std::cerr << "Unexpected close of sender" << err;
          delete[] buff;
          return false;
      }
      totalRecvd += static_cast<size_t>(recvd);

      // send received bytes to receiver
      sendSuccess = receiver->sendall(buff, static_cast<size_t>(recvd));
      if (!sendSuccess) {
        std::cerr << "Unexpected fail of receiver";
        delete[] buff;
        return false;
      }
    }
  } else { // tunnel until sender disconnects
    while (recvd > 0 && sendSuccess)
    {
      // receive a max of chunksize bytes
      recvd = recv(sendFd, buff, chunksize, 0);
      err = errno;
      if (recvd == -1) {
          std::cerr << "Unexpected fail of sender" << err;
          delete[] buff;
          return false;
      }

      // send received bytes to receiver
      sendSuccess = receiver->sendall(buff, static_cast<size_t>(recvd));
      if (!sendSuccess) {
        std::cerr << "Unexpected fail of receiver";
        delete[] buff;
        return false;
      }
    }
  }
  std::cout << std::endl;
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
