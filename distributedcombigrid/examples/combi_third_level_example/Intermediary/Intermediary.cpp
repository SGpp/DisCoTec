#include "Communicator.hpp"
#include "Intermediary.hpp"
#include "../../../src/sgpp/distributedcombigrid/third_level/NetworkUtils.hpp"

Intermediary::Intermediary(int port) : port_(port) {
}

Intermediary::~Intermediary() {
  for (int i = 0; i < comms_.size(); i++) {
    delete comms_[i];
  }
}

void Intermediary::runForever() {
  acceptClient();
}

/*
 * Continually listens on port for new clients. If a client connects, it
 * either gets accepted and becomes a new participant of a communicator or an
 * already existing participant wants to create a new connection i.e. for
 * data transfer.
 */
void Intermediary::acceptClient() {
  ServerSocket server(port_);
  for (;;)
  {
    ClientSocket* client = server.acceptClient();
    std::thread handle(&Intermediary::handleClient, this, client);
    handle.detach();
  }
}

void Intermediary::handleClient(ClientSocket* client) {
  std::vector<std::string> argv;
  std::string initmsg;
  assert(client->isInitialized() && "Accepted client corrupt");
  client->recvallPrefixed(initmsg);
  NetworkUtils::split(initmsg, '#', argv);

  if (argv[0] == "p") { // new participant.
    assert(argv.size() == 2 &&
        "New participants first message must be of form p#commRank");
    assert(NetworkUtils::isInteger(argv[1]) &&
        "New participants desired comm rank must be a number");
    int commID = std::stoi(argv[1]);

    std::cout << "Accepted new client with init mesg: " << initmsg << "\n";

    if (comms_.find(commID) != comms_.end()) {
      Communicator* comm = comms_[commID];
      comm->addParticipant(client);
    } else {
      Communicator* comm = new Communicator(*this, commID);
      comms_.insert(std::pair<int, Communicator*>(commID, comm));
      comm->addParticipant(client);
    }
  }
  else if (argv[0] == "d") { // new data connection.
    std::cout << "new data request: " << initmsg << std::endl;
    assert(argv.size() == 3 &&
        "Request for new data connection must be of form d#commRank#partRank");
    assert(NetworkUtils::isInteger(argv[1]) &&
      "Participants comm id must be a number");
    int commID = std::stoi(argv[1]);

    assert(NetworkUtils::isInteger(argv[2]) &&
      "Participants rank must be a number");
    int rank = std::stoi(argv[2]);

    assert(comms_.find(commID) != comms_.end());
    Communicator* comm = comms_[commID];
    assert(comm->existsParticipant(rank));

    Participant& part = comm->getParticipant(rank);
    {
      std::lock_guard<std::mutex> lock(part.newConnLock_);
      part.newConn_ = client;
      part.newConnValid_ = true;
      part.newConnCV_.notify_one();
    }
  } else { // unknown participant
    assert(false && "Initial message not valid");
    delete client;
  }
}
