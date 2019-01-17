#ifndef MESSAGEUTILS_HPP
#define MESSAGEUTILS_HPP

#include "SimpleAmqpClient/SimpleAmqpClient.h"

const int   noTimeout = -1;

class MessageUtils
{
  public:
    static void sendMessage(const std::string&               message,
                            const std::string                queue,
                            const AmqpClient::Channel::ptr_t channel);

    static bool receiveMessage(const AmqpClient::Channel::ptr_t channel,
                               const std::string                queue,
                               std::string&                     message,
                               int                              timeout);
};

#endif
