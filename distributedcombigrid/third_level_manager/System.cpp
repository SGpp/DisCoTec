#include "System.hpp"

namespace combigrid {

System::System(std::shared_ptr<ClientSocket>& connection, size_t id)
  : connection_(connection), id_(id)
{
}

void System::sendMessage(const std::string& message) const
{
  connection_->sendallPrefixed(message);
}

bool System::receiveMessage(std::string& message, int timeout) const
{
  return connection_->recvallPrefixed(message);
}

bool System::receivePosNumber(size_t& number) const
{
  std::string message;
  if (receiveMessage(message))
  {
    assert(NetworkUtils::isInteger(message));
    number = std::stoul(message);
    return true;
  }
  return false;
}

bool System::hasMessage(int timeout) {
  return connection_->isReadable(timeout);
}

std::shared_ptr<ClientSocket> System::getConnection() const
{
  return connection_;
}

size_t System::getId() const
{
  return id_;
}

}
