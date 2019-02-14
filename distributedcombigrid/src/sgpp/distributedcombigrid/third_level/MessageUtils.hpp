#ifndef MESSAGEUTILS_HPP
#define MESSAGEUTILS_HPP

#include "SimpleAmqpClient/SimpleAmqpClient.h"


class MessageUtils
{
  public:
    static const int noTimeout = -1;


    static std::string createMessageQueue(const std::string&               name,
                                          const AmqpClient::Channel::ptr_t channel);

    static void sendMessage(const std::string&               message,
                            const std::string                queue,
                            const AmqpClient::Channel::ptr_t channel);

    static bool receiveMessage(const AmqpClient::Channel::ptr_t channel,
                               const std::string                consumerTag,
                               const std::string                queue,
                               std::string&                     message,
                               int                              timeout = noTimeout);

    static std::string setupConsumer(const AmqpClient::Channel::ptr_t channel,
                                     const std::string                queue);
};

#endif
