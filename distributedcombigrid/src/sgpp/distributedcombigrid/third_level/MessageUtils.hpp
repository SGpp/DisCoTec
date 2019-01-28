#ifndef MESSAGEUTILS_HPP
#define MESSAGEUTILS_HPP

#include "SimpleAmqpClient/SimpleAmqpClient.h"


class MessageUtils
{
  private:
    static const int noTimeout = -1;

  public:
    static void sendMessage(const std::string&               message,
                            const std::string                queue,
                            const AmqpClient::Channel::ptr_t channel);

    static bool receiveMessage(const AmqpClient::Channel::ptr_t channel,
                               const std::string                queue,
                               std::string&                     message,
                               int                              timeout = noTimeout);
};

#endif
