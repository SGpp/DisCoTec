#include "System.hpp"

System::System(std::shared_ptr<ClientSocket>& connection)
  : connection_(connection)
{
}

void System::sendMessage(const std::string& message) const
{
  connection_->sendallPrefixed(message);
}

bool System::receiveMessage(std::string& message, int timeout) const
{
  if (connection_->isReadable(timeout))
  {
    connection_->recvallPrefixed(message);
    return true;
  }
  return false;
}

bool System::receivePosNumber(size_t& number, int timeout) const
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

std::shared_ptr<ClientSocket> System::getConnection() const
{
  return connection_;
}
