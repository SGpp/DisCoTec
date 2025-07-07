#include "discotec/third_level/NetworkUtils.hpp"
#include <string>
#include <thread>

namespace combigrid {

class System
{
  private:
    std::shared_ptr<ClientSocket> connection_;
    size_t id_;

  public:
    System(std::shared_ptr<ClientSocket>& connection, size_t id);

    void sendMessage(const std::string& message) const;
    bool receiveMessage(std::string& message, int timeout=NetworkUtils::noTimeout) const;
    bool receivePosNumber(size_t& number) const;
    bool hasMessage(int timeout);
    size_t getId() const;

    std::shared_ptr<ClientSocket> getConnection() const;
};

}
