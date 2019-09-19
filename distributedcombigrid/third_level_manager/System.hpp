#include "sgpp/distributedcombigrid/third_level/NetworkUtils.hpp"
#include <string>
#include <thread>

class System
{
  private:
    std::shared_ptr<ClientSocket> connection_;

  public:
    System(std::shared_ptr<ClientSocket>& connection);

    void sendMessage(const std::string& message) const;
    bool receiveMessage(std::string& message, int timeout=NetworkUtils::noTimeout) const;
    bool receivePosNumber(size_t& number, int timeout=NetworkUtils::noTimeout) const;

    std::shared_ptr<ClientSocket> getConnection() const;
};
